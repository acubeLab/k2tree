#!/usr/bin/env python3

import sys, argparse

Description = """
Convert a full binary matrix in ascii to sparse format by considering 
the 0/1 characters as matrix entries and ignoring all other characters (e.g., spaces, new lines).
Each line of the input file is processed sequentially and represent a matrix row
"""

def main():
  parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('input', help='input matrix file name', type=str)
  parser.add_argument("-r", help="reverse conversion", action="store_true")
  parser.add_argument("-t", help="transpose matrix", action="store_true")
  parser.add_argument('-o', metavar='outfile', help='output file name (def. input.sparse/output.full)',type=str,default="" )
  args = parser.parse_args()

  if args.r:
    # reverse conversion: from sparse to full binary matrix
    totnz, rows,cols = reverse_conversion(args)
  else:
    # convert from full binary matrix to sparse format
    with open(args.input,"rt") as f:
      # set output file name and output file mode
      if args.o!="": 
        outname = args.o
      else:
        outname = f"{args.input}.sparse"
      totnz = 0
      rindex = 0
      maxcindex = -1
      with open(outname,"w") as g:
        for line in f:
          cindex = 0
          for c in line:
            if c=='0' or c=='1':
              if c=='1':
                if args.t:
                  g.write(f"{cindex} {rindex}\n") # transpose
                else:
                  g.write(f"{rindex} {cindex}\n")
                totnz += 1
                maxcindex = max(maxcindex,cindex)
              cindex += 1
          rindex += 1
    rows = rindex
    cols = maxcindex+1
  print("Number of matrix entries written:", totnz)
  print("Matrix dimensions: ", rows, "x", cols)
  print("==== Done")


# reverse conversion: from sparse to full binary matrix
# not very efficient 
def reverse_conversion(args):
    # set output file name and output file mode
    if args.o!="":  outname = args.o
    else:           outname = args.input + ".full"
    # first pass: determine matrix dimensions
    maxrindex = -1
    maxcindex = -1
    nonzeros = {}
    with open(args.input,"rt") as f:
      # read nonzeros
      for line in f:
        items = line.split()
        if len(items)!=2:
          continue
        rindex = int(items[0])
        cindex = int(items[1])
        if args.t:
          rindex,cindex = cindex,rindex  # transpose
        maxrindex = max(maxrindex, rindex)
        maxcindex = max(maxcindex, cindex)
        if rindex not in nonzeros:
          nonzeros[rindex] = set()
        nonzeros[rindex].add(cindex)
    # second pass: write full matrix
    totnz = 0
    with open(outname,"wt") as g:
      for r in range(maxrindex+1):
        row = ""
        for c in range(maxcindex+1):
          if r in nonzeros and c in nonzeros[r]:
            row += '1'
            totnz += 1
          else:
            row += '0'
        g.write(row +'\n')
    return totnz,maxrindex+1,maxcindex+1

if __name__ == '__main__':
    main()
