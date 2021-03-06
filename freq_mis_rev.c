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

char g_idx_dict['T' + 1];
typedef enum Base { A, C, G, T, KMER_SIZE = 2 } Base;
void init_dict()
{
  g_idx_dict['A'] = A;
  g_idx_dict['C'] = C;
  g_idx_dict['G'] = G;
  g_idx_dict['T'] = T;
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

static int g_max;
void inc_neighbours(unsigned int kmer, unsigned int k, unsigned int d, unsigned int *tmap)
{
  ++tmap[kmer];
  if (tmap[kmer] > g_max)
    g_max = tmap[kmer];
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
	      {
		++tmap[flip];
		if (tmap[flip] > g_max)
		  g_max = tmap[flip];
	      }
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
  int k, d;
  sscanf(vars, "%d %d", &k, &d);
  --vars;
  char *genome = buf;


  size_t tmap_size = pow(4, k) * sizeof(int);

  int *tmap = (int *)malloc(tmap_size);
  memset(tmap, 0, tmap_size);

  printf("k = %d, d = %d\n", k, d);
  init_dict();
  g_max = 0;

  int idx = 0;
  char *cp = buf;
  int i;

  for (i = 0; i < k; ++i)
    {
      idx = (idx << 2) + g_idx_dict[*cp];
      ++cp;
    }

  int kmask = (1 << (2*k)) - 1;
  inc_neighbours(idx, k, d, tmap);
  int reverse_idx = (~idx) & kmask;
  inc_neighbours(reverse_idx, k, d, tmap);

  while (*cp != '\n')
    {
      idx = ((idx << 2) & kmask) + g_idx_dict[*cp];
      inc_neighbours(idx, k, d, tmap);

      reverse_idx = (~idx) & kmask;
      inc_neighbours(reverse_idx, k, d, tmap);
      //      inc_
      ++cp;
    }


  for (i = 0; i < tmap_size; ++i)
    {
      if (tmap[i] == g_max)
	{
	  print_clump(i, k);
	}
    }

  free(tmap);
  return 0;
}
