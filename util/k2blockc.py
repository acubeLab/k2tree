#!/usr/bin/env python3

import sys, argparse, subprocess, os, concurrent.futures

K2TOOL = "../k2sparse.x"
K2EXT = ".k2"

Description = f"""
Takes as input a text file in .rowm format containing pairs of row/column 
indices of the nonzero elements sorted by row index and a parameter b.
Outputs b files in the same textual format each containing
a set of consecutive rows. The splitting is done dividing almost 
equally the number of nonzeros among the blocks.
The files are later compressed individually invoking {K2TOOL}

The nonzeros are almost equally distributed since the nonzers of a given row 
always go together in the same file. Note that the index of the elements 
are not modified in the copy so formally the output matrices all have the same
size of the input matrix even if most of the rows will be empty. 

If the number of nonzero is not provided with the -n option the program will
get it counting the number of lines in the input file.

The input file is assumed to be simply sorted by row: it is not necessary that 
elements in the same row are also sorted by column"""



def main():
  parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('input', help='input matrix file name', type=str)
  parser.add_argument('blocks', help='number of input blocks', type=int)
  parser.add_argument('-n', help='number of nonzeros (def. number of input file lines)',type=int,default=-1)
  parser.add_argument('-s', help='matrix size (def. # of int32 in .ccount file)',type=int,default=-1)
  parser.add_argument('-o', metavar='outfile', help='output file base name (def. input file name)',type=str,default="" )
  parser.add_argument('-k', help="keep temporary files", action='store_true')
  parser.add_argument('--cext', help=f"extension for compressed files (def. {K2EXT})", default=K2EXT,type=str)
  parser.add_argument('--copts', help=f"compression options for {K2TOOL} (def. none)", default='',type=str)
  parser.add_argument('-v', help="verbose", action='store_true')
  args = parser.parse_args()
    
  if args.blocks<1:
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
  # get matrix size if missing
  if args.s < 0:
    ccount_file = os.path.splitext(args.input)[0]+'.ccount'
    if not os.path.exists(ccount_file):
      print(f"Error: no matrix size provided and file {ccount_file} not found")
      sys.exit(1)
    args.s = os.path.getsize(ccount_file)//4
  assert args.s>0, f"Error: invalid matrix size: {args.s}" 
  if args.v:
    print("Creating matrices of size:",args.s)
  with open(args.input,"rt") as f:  
    # get number of nonzeros from the input file if missing 
    if args.n<0: 
      args.n = sum(1 for _ in f)
      f.seek(0) # rewind input file
    if args.v:
      print("Number of nonzeros in the input file", args.n) 
    # split the input file into blocks
    tot_lines = 0
    tot_written = 0
    lastrow = -1
    lastline = ""
    with concurrent.futures.ThreadPoolExecutor() as executor:
      futures = []
      for i in range(args.blocks):
        # create and compress i-th text file  
        gname = f"{args.input}.tmp.{args.blocks}.{i}"
        outname = f"{args.o}.{args.blocks}.{i}{args.cext}"
        with open(gname,"wt") as g:
          target = (args.n-tot_written)//(args.blocks-i)
          assert target>0, f"Error: not enough nonzeros to split in {args.blocks} blocks"  
          if args.v: print(f"  Aiming at {target} nonzeros in file {gname}")
          # write eventual suspended line from previous block
          if len(lastline)>0:
            g.write(lastline)  # write pending line
            tot_written += 1
            target -= 1
          # copy target lines from f to g
          for _ in range(target):
            line = f.readline()
            tot_lines += 1
            row = int(line.split()[0])
            assert(row>=lastrow), f"Error: input file not sorted by row at line {tot_lines}: {line}"
            lastrow = row
            g.write(line)
            tot_written += 1
          # complete current row
          while True:
            line = f.readline()
            if len(line)==0: break
            tot_lines += 1
            row = int(line.split()[0])
            assert(row>=lastrow), f"Error: input file not sorted by row at line {tot_lines}: {line}"
            if row>lastrow: 
              lastrow = row
              lastline = line
              break
            g.write(line)
            tot_written += 1
        # gname completed: compress and delete it 
        futures.append(executor.submit(k2compress,gname,outname,args))
    # wait for all futures to complete and report errors  
    for f in concurrent.futures.as_completed(futures):
      if f.result(): print(f.result())  
    # pool operations completed, check if all nonzeros were written    
  assert tot_lines==args.n, f"Mismatch in nonzeros: expected {args.n}, written {tot_lines}" 
  assert tot_written==args.n, f"Mismatch in nonzeros: expected {args.n}, written {tot_written}" 
  print("==== Done")



# compress a temporary file and then if requested delete it
def k2compress(fname,outname,args):
  cmd = f"{args.exe} -o {outname} -s {args.s} {args.copts} {fname}"
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
