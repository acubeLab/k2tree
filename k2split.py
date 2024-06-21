#!/usr/bin/env python3

import sys, argparse, subprocess, os, concurrent.futures

K2TOOL = "k2sparse.x"

Description = """
Quick and dirty tool to take a k2 compressed file and split it into an assigned number of
smaller k2 compressed files containing roughly the same number of nonzeros.

Relies of the textual output of k2sparse.x to determine the number of nonzeros in the input file."""


def main():
  parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('input', help='input matrix file name', type=str)
  parser.add_argument('blocks', help='number of input blocks', type=int)
  parser.add_argument('-o', metavar='outfile', help='output file base name (def. input file name)',type=str,default="" )
  parser.add_argument('-k', help="keep temporary files", action='store_true')
  parser.add_argument('-opts', help=f"compression options for {K2TOOL} (def none)", default='',type=str)
  parser.add_argument('-v', help="verbose", action='store_true')
  args = parser.parse_args()
    
  if args.blocks<=1:
    print("The number of blocks must be at least 1")
    sys.exit(1)
  #  base for output file name   
  if args.o=="":  args.o = args.input
  # temporary decompressed file   
  tmpname = args.o+".tmpfile"
  # get path to k2split.py directory 
  args.main_dir = os.path.dirname(sys.argv[0])
  args.exe = os.path.join(args.main_dir,K2TOOL)

  # decompress k2 file and recover # nonzeros from stdout
  result = subprocess.run(f"{args.exe} {args.input} -dv -o {tmpname}".split(), 
           capture_output=True, text=True)
  if result.returncode!=0:
    print("Decompression failed with exit code:", result.returncode)
    print("stderr from the decompression program:")
    print(result.stderr)
    print("==== Exiting")
    sys.exit(1)
  # decompression ok, retrieve # of nonzeros from stdout                    
  data = result.stdout.split("\n")[2]
  nonz = int(data.split(",")[4].split()[-1]) 
  if args.v:
    print("stdout from the decompression program:")
    print(result.stdout)     
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
        futures.append(executor.submit(k2compress,gname,args,i))
      # wait for all futures to complete and report errors  
      for f in concurrent.futures.as_completed(futures):
        if f.result(): print(f.result())  
    # pool operations completed, check if all nonzeros were written    
  assert check_lines==nonz, f"Mismatch in nonzeros: expected {nonz}, written {check_lines}" 
  print("==== Done")

# compress a temporary file and then if requested delete it
def k2compress(fname,args,i):

  cmd = f"{args.exe} -o {args.o}.{args.blocks}.{i} {args.opts} {fname}"
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
