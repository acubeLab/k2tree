#!/usr/bin/env python3

import sys, argparse, subprocess, os, concurrent.futures

K2TOOL = "k2sparse.x"

Description = """
Quick and dirty tool that takes a text file containing the list of pairs <row col> 
of the positions of the nonzeros of a binary matrix in row major order, 
and a parameter b, and splits the list into b files (almost) equally distributing 
the nonzeros among them. The files are later compressed individually with k2sparse.x

The nonzeros are almost equally distributed since the nonzers of a given row always 
go together in the same file: we are splitting the matrix by rows since this 
simplify the matrix-vector multiplication. Note that the index of the row elements 
are not modified in the copy.

If the number of nonzero is not provided with the -n option the program will
get it counting the number of lines in the input file.

The input file is assumed to be sorted by row: it is not necessary that elements
in the same row are also sorted by column"""



def main():
  parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('input', help='input matrix file name', type=str)
  parser.add_argument('blocks', help='number of input blocks', type=int)
  parser.add_argument('-N', help='number of nonzeros (def. number of input file lines)',type=int,default=-1)
  parser.add_argument('-o', metavar='outfile', help='output file base name (def. input file name)',type=str,default="" )
  parser.add_argument('-k', help="keep temporary files", action='store_true')
  parser.add_argument('--copts', help=f"compression options for {K2TOOL} (def none)", default='',type=str)
  parser.add_argument('-v', help="verbose", action='store_true')
  args = parser.parse_args()
    
  if args.blocks<=1:
    print("The number of blocks must be at least 1")
    sys.exit(1)
  #  base for output file name   
  if args.o=="":  args.o = args.input
  # get path to k2split.py directory 
  args.main_dir = os.path.dirname(sys.argv[0])
  args.exe = os.path.join(args.main_dir,K2TOOL)
  if not os.path.exists(args.exe):
    print(f"Error: {args.exe} not found")
    sys.exit(1)
  # split the input file into blocks  
  with open(args.input,"rt") as f:  
    # get number of nonzeros from the input file if missing 
    if args.n<0: 
      args.n = sum(1 for _ in f)
      f.seek(0) # rewind input file
    if args.v:
      print("Number of nonzeros in the input file", args.n) 
    # split the input file into blocks
    remaining = args.n
    tot_lines = 0
    lastrow = -1
    lastline = ""
    for i in range(args.blocks):
      # create and compress i-th text file  
      gname = f"{args.input}.tmp.{args.blocks}.{i}"
      outname = f"{args.o}.{args.blocks}.{i}"
      with open(gname,"wt") as g:
        target = remaining//(args.blocks-i)
        if args.v:
          print(f"  Aiming at {target} nonzeros in file {gname}")
        # copy target lines from f to g
        if len(lastline)>0:
          g.write(lastline)  # write pending line
          target -= 1
        for _ in range(target):
          line = f.readline()
          tot_lines += 1
          row = int(line.split()[0])
          assert(row>=lastrow), f"Error: input file not sorted by row at line {tot_lines}: {line}"
          lastrow = row
          g.write(line)
      # compress gname and delete it 
      k2compress(gname,outname,args)
    



  # decompress k2 file and recover # nonzeros from stdout
  result = subprocess.run(f"{args.exe} {args.input} -dv -o {tmpname}".split(), 
           capture_output=True, text=True)
  if result.returncode!=0:
    print("Decompression failed with exit code:", result.returncode)
    print("stderr from the decompression program:")
    print(result.stderr)
    print("==== Exiting")
    sys.exit(1)
  # decompression ok, retrieve size and # of nonzeros from stdout 
  # print(result.stdout) 
  lines = result.stdout.split("\n")  #  lines of stdout
  args.size = int(lines[2].split(",")[0].split()[-1]) # get compressed matrix size              
  nonz = int(lines[3].split(",")[4].split()[-1])     # get number of nonzeros
  if args.v:
    print("stdout from the decompression program:")
    print(result.stdout) 
  print("Size of the input matrix", args.size)        
  print("Number of nonzeros in the input file", nonz)
  if nonz<args.blocks:
    print("Nothing to do: Number of nonzeros smaller than the number of blocks!")
    print("==== Exiting")
    if not args.k: os.unlink(tmpname)
    sys.exit(1)
  # copy tmpfile content to blocks distinct files  
  with open(tmpname,"rt") as f:
    if not args.k: os.unlink(tmpname)   #  no longer needed
    check_lines=0
    with concurrent.futures.ThreadPoolExecutor() as executor:
      futures = []
      for i in range(args.blocks):
        # create i-th text file 
        gname = f"{tmpname}.{args.blocks}.{i}"
        outname = f"{args.o}.{args.blocks}.{i}"
        with open(gname,"wt") as g:
          totlines = nonz//args.blocks
          if i < nonz%args.blocks: totlines += 1
          check_lines += totlines
          if args.v:
            print(f"  Writing {totlines} nonzeros to file {gname}") 
          for _ in range(totlines):
            r = f.readline()
            g.write(r)
        # compress gname and delete it 
        futures.append(executor.submit(k2compress,gname,outname,args))
      # wait for all futures to complete and report errors  
      for f in concurrent.futures.as_completed(futures):
        if f.result(): print(f.result())  
    # pool operations completed, check if all nonzeros were written    
  assert check_lines==nonz, f"Mismatch in nonzeros: expected {nonz}, written {check_lines}" 
  print("==== Done")

# compress a temporary file and then if requested delete it
def k2compress(fname,outname,args):

  cmd = f"{args.exe} -o {outname} -s {args.size} {args.copts} {fname}"
  if args.v:
    print("  Executing:",cmd)
  result = subprocess.run(cmd.split(), 
           capture_output=True, text=True)
  if result.returncode!=0:
    print("  Compression failed with exit code:", result.returncode)
    print("  stderr from the compression program:")
    print("  "+result.stderr)
  elif args.v and len(result.stdout)>0:
    print("  stdout from the decompression program:")
    print("  "+result.stdout)     
  if not args.k: os.unlink(fname)


if __name__ == '__main__':
    main()
