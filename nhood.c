#include <stdio.h>
#include <fcntl.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


int open_file(int argc, char *argv[])
{
  if (argc < 2)
    {
      printf("Usage: %s filename\n", argv[0]);
      return 0;
    }
  int fd = open(argv[1], O_RDONLY);
  if (fd == -1)
    printf("Can't open file %s\n", argv[1]);
  return fd;
}

int read_file(int fd, char *buf, int buf_size)
{
  if (fd == -1)
    return 0;

  int len = read(fd, buf, buf_size);

  if (len == -1)
    {
      printf("Error reading file\n");
      exit(0);
    }
  if (len == buf_size)
    {
      printf("File is too big\n");
      exit(0);
    }
  return len;
}

//returns line start
char * find_line_start(char *buf_start, char *buf_end)
{
  while (*buf_end != '\n')
    {
      --buf_end;
      if (buf_start == buf_end)
	{
	  printf("Wrong file format\n");
	  exit(0);
	}
    }
  ++buf_end;
  return buf_end;
}

enum { BUF_SIZE = 5*1024*1024 }; //x mb


// 2 bits representation
typedef enum Base { A, C, G, T, KMER_SIZE = 2 } Base;
Base calc(char ch)
{
  switch (ch)
    {
    case 'A':
      return A;
    case 'C':
      return C;
    case 'G':
      return G;
    case 'T':
      return T;
    }
  printf("Wrong input %c\n", ch);
  exit(0);
}

void print_clump(int kmer_idx, int k)
{
  char buf[k+1];
  buf[k] = 0;
  static const char labels[] = { 'A', 'C', 'G', 'T' };
  while (k--)
    {
      buf[k] = labels[kmer_idx & 0x03];
      kmer_idx >>= 2;
    }
  printf("%s\n", buf);
}


int get_flip(unsigned int flip, unsigned int dmask, unsigned int flipmap)
{
  int flip_out = flip;
  //  printf("dmask %x, flipmap %x, ", dmask, flipmap);
  unsigned int flipbitpos = 0;
  unsigned int dbit = (dmask ^ (dmask - 1));
  dbit = dbit ^ (dbit >> 1); // get first one-bit
  while (dbit)
    {
      //      getchar();
      const unsigned int shift = dbit * dbit;
      //      printf("dbit %x, shift %x \n", dbit, shift);
      const unsigned int kmerc = flip & (0x03 * shift);
      const unsigned int flipshift = flipbitpos*2;
      const unsigned int flipbit = 0x03 << flipshift;
      const unsigned int flipc = ((flipmap & flipbit) >> flipshift) * shift;

      /* print_clump(kmerc, 3); */
      /* print_clump(flipc, 3); */
      /* getchar(); */

      if (kmerc == flipc) //no flip this flipmap, give me next
	return flip_out;
      flip &= ~(0x03 * shift);
      flip |= flipc;

      dmask = dmask & ~dbit;
      dbit = dmask;
      if (dbit == 0)
	break;
      dbit = (dbit ^ (dbit - 1)); //  and find next
      dbit = dbit ^ (dbit >> 1);
      //      printf("dbit %x get next\n", dbit);
      ++flipbitpos;
    }
  return flip;
}

void inc_neighbours(unsigned int kmer, unsigned int k, unsigned int d, unsigned int *tmap)
{
  print_clump(kmer, k);
  unsigned int dstop = 1 << k;
  while (d > 0)
    {
      unsigned int dmask = (1 << d) - 1;
      unsigned int flip_stop = (1 << 2*d);
      while (dmask < dstop)
	{
	  unsigned int flip_map = 0;
	  for (;flip_map < flip_stop; ++flip_map)
	  {
	    unsigned int flip = get_flip(kmer, dmask, flip_map);
	    if (flip != kmer)
	      print_clump(flip, k);
	  }
	  unsigned int t = dmask | (dmask - 1); 
	  // t gets v's least significant 0 bits set to 1
	  // Next set to 1 the most significant bit to change, 
	  // set to 0 the least significant ones, and add the necessary 1 bits.
	  dmask = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(dmask) + 1)); 
	}
      --d;
    }
}

int main(int argc, char **argv)
{
  int fd = open_file(argc, argv);

  char buf[BUF_SIZE];
  int len = read_file(fd, buf, sizeof(buf));

  char *vars = find_line_start(buf, buf + len - 2);
  int d;
  sscanf(vars, "%d", &d);

  int idx = 0;
  char *cp = buf;
  int k = 0;
  while (*cp != '\n')
    {
      idx = (idx << 2) + calc(*cp);
      ++cp;
      ++k;
    }
  printf("k = %d, d = %d\n", k, d);

  int k_mask = (1 << (2*k)) - 1;
  inc_neighbours(idx, k, d, NULL);

  return 0;
}
