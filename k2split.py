#!/usr/bin/env python3

import sys, argparse, subprocess, os

K2TOOL = "./k2sparse.x"



Description = """
Delete from a textual matrix in one-entry-per-line format
all entries involving an index with id >= newsize  
"""

def main():
  parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('input', help='input matrix file name', type=str)
  parser.add_argument('blocks', help='number of input blocks', type=int)
  parser.add_argument('-o', metavar='outfile', help='output file base name (def. input file name)',type=str,default="" )
  parser.add_argument('-v', help="verbose", action='store_true')
  args = parser.parse_args()
    
  if args.blocks<=1:
    print("The number of blocks must be at least 1")
    sys.exit(1)
  #  base for output file name   
  if args.o!="": 
    outname = args.o
  else:
    outname = args.input
  # temporary decompressed file   
  tmpname = outname+".tmpfile"

  # decompress k2 file and recover # nonzeros from stdout
  result = subprocess.run(f"{K2TOOL} {args.input} -dv -o {tmpname}".split(), 
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
    os.unlink(tmpname)
    sys.exit(1)
  # copy tmpfile content to blocks distinct files  
  with open(tmpname,"rt") as f:
    os.unlink(tmpname)   #  nolonger needed
    check_lines=0
    for i in range(args.blocks):
      outi = f"{outname}.{args.blocks}.{i}"
      with open(outi,"wt") as g:
        totlines = nonz//args.blocks
        if i < nonz%args.blocks: totlines += 1
        check_lines += totlines
        if args.v:
          print(f"  Writing {totlines} nonzeros to file {outi}") 
        for _ in range(totlines):
          r = f.readline()
          g.write(r)
  print("Total nonzero written:", check_lines)
  assert check_lines==nonz, f"Mismatch in nonzeros: expected {nonz}, written {check_lines}" 
  # compress partial files 

  
  print("==== Done")


if __name__ == '__main__':
    main()
