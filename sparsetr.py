#!/usr/bin/env python3

import sys, argparse

Description = """
Transform a sparse matrix (one-line per entry format) according to the options:
 - force matrix size (discarding elements with index >= size)
 - transpose
 - set to 1 the diagonal entries
 - add elements to make the matrix symmetric
 The output is again in sparse format (one-line per entry)
"""

def main():
  parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('input', help='input matrix file name', type=str)
  parser.add_argument("-s", metavar="size", help="force matrix size", type=int,default=-1)
  parser.add_argument("-t", help="transpose matrix", action="store_true")
  parser.add_argument("-d", help="set to 1 the diagonal entries", action="store_true")
  parser.add_argument("-S", help="make matrix symmetric", action="store_true")
  parser.add_argument('-o', metavar='outfile', help='output file name (def. input.trx)',type=str,default="" )
  args = parser.parse_args()

  # set output file name and output file mode
  if args.o!="":  outname = args.o
  else:           outname = args.input + ".trx"

  if args.t  and args.S:
    print("Options -t and -S together do not make sense!",file=sys.stderr)
    exit(1)

  maxrindex = -1
  maxcindex = -1
  diagonal = set()
  offdiagonal = set()  # used to make symmetric
  totnz = 0
  with open(outname,"wt") as g:
    with open(args.input,"rt") as f:
      # read nonzeros
      for line in f:
        if line[0]=='#':  # skip comment lines
          continue
        items = line.split()
        if len(items)!=2:   # skip malformed lines
          continue
        rindex = int(items[0])
        cindex = int(items[1])
        if args.t:
          rindex,cindex = cindex,rindex  # transpose
        if args.s>=0:
          if rindex>=args.s or cindex>=args.s:
            continue  # skip entries outside the forced size
        maxrindex = max(maxrindex, rindex)
        maxcindex = max(maxcindex, cindex)
        if cindex == rindex:
          diagonal.add(rindex)  # keep track of existing diagonal entries
        elif args.S:
          offdiagonal.add( (rindex,cindex) )
        totnz += 1
        g.write(f"{rindex} {cindex}\n")
    # add diagonal entries if needed
    if args.d:
      size = args.s if args.s>=0 else max(maxrindex,maxcindex)+1
      for i in range(size):
        if i not in diagonal:
          totnz += 1
          g.write(f"{i} {i}\n")
    # add offdiagonal to make symmetric
    if args.S:
      for (r,c) in offdiagonal:
        if (c,r) not in offdiagonal:
          # add element (c,r)
          maxrindex = max(maxrindex, c)
          maxcindex = max(maxcindex, r)
          totnz += 1
          g.write(f"{c} {r}\n")
  print("Number of matrix entries written:", totnz)
  print("Actual matrix dimensions: ", maxrindex+1, "x", maxcindex+1)
  return

if __name__ == '__main__':
    main()
