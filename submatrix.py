#!/usr/bin/env python3

import sys, argparse

Description = """
Delete from a textual matrix in one-entry-per-line format
all entries involving an index with id >= newdime  
"""

def main():
  parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('input', help='input matrix file name', type=str)
  parser.add_argument('newsize', help='output matrix size', type=int)
  parser.add_argument('-o', metavar='outfile', help='output file name (def. input.newsize)',type=str,default="" )
  args = parser.parse_args()
    
  if args.newdim<=1:
    print("The new matrix size must be at least 1")
    sys.exit(1)   
    
  with open(args.input,"rt") as f:
    # set output file name and output file mode
    if args.o!="": 
      outname = args.o
    else:
      outname = f"{args.input}.{args.newdim}"
    maxinput = maxoutput = -1  
    with open(outname,"w") as g:
      r = 0 # number of read rows
      wr = 0 # number of written rows
      for line in f:
        r += 1
        p = line.split()
        if len(p)!=2:
          print("Illegal entry at line", r)
          sys.exit(1)
        a = int(p[0]); b = int(p[1])
        if maxinput < 0: maxinput = max(a,b)
        else: maxinput = max(maxinput,a,b)
        if a<args.newdim and b < args.newdim:
          if maxoutput < 0: maxoutput = max(a,b)
          else: maxoutput = max(maxoutput,a,b)
          print(f"{a} {b}",file=g)
          wr +=1
  print("Number of input entries:", r)
  print("Input matrix size:", maxinput, "Density:", r/(maxinput*maxinput))
  print("Number of output entries:", wr)
  print("Output matrix size:", maxoutput, "Density:", wr/(maxoutput*maxoutput))
  print("==== Done")


if __name__ == '__main__':
    main()
