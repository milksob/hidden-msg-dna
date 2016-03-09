#include <stdio.h>
#include <fcntl.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

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

static int gmark;
void mark_neighbours(unsigned int kmer, unsigned int k, unsigned int d, unsigned int *tmap)
{
  //  print_clump(kmer, k);
  tmap[kmer] |= gmark;
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
		tmap[flip] |= gmark;
		//       		print_clump(flip, k);
		//		getchar();
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

int calc_distance(int kmer_idx, int pattern_idx, int k)
{
  //    print_clump(kmer_idx, k);
  //    print_clump(pattern_idx, k);

  static int diff;
  diff = kmer_idx ^ pattern_idx;
  static const int odd_mask = 0xAAAAAAAA;   //10 10 ... mask
  static const int even_mask = 0x55555555;  //01 01 ... mask

  diff |= ((diff & odd_mask) >> 1); // OR evens with odds
  diff &= even_mask; // clear all odds;
  diff = __builtin_popcount(diff);  // number of different letters

  //    printf("calc_distance: %d\n", diff);
    //    getchar();
  return diff;
}

float calc_prob(int idx, int k, float **prf)
{
  float prob = 1;
  int i, c;
  //  print_clump(idx, k);
  for (i = 0; i < k; ++i)
    {
      c = 2*(k - 1 - i);
      c = idx >> c;
      c &= 0x03;
      //      const char sym[] = {'a','c','g','t'};
      prob = prob * prf[c][i];
      //      printf(".. %c %d prf = %f\n", sym[c], i, prf[c][i]);
      //            getchar();
    }
  return prob;
}

int main(int argc, char **argv)
{
  int fd = open_file(argc, argv);

  char buf[BUF_SIZE];
  int len = read_file(fd, buf, sizeof(buf));

  float *prf[4]; // symbols a c g t
  int k, i, j;
  char *cp = buf;// + (len - 4*16 - 3 - 10);
  while (*cp != '\n')
    ++cp;
  ++cp;

  k = strtol(cp, &cp, 10);
  //    printf("k = %d\n", k);
  for (i = 0; i < 4; ++i)
    {
      prf[i] = (float*)malloc(k*sizeof(float));
      for (j = 0; j < k; ++j)
	{
	  prf[i][j] = strtod(cp, &cp);
	  //	  printf("prf[%d][%d] <- %f\n", i, j, prf[i][j]);
	}
    }

  init_dict();

  size_t tmap_size = pow(4, k) * sizeof(int);
  int kmask = (1 << (2*k)) - 1;

  int idx = 0;
  cp = buf;
  for (j = 0; j < k; ++j)
    {
      if (*cp == '\n')
	assert(0 && "wrong format file");
      idx = (idx << 2) + g_idx_dict[*cp];
      ++cp;
    }
  float prob = calc_prob(idx, k, prf);
  float prob_max = prob;
  float prob_max_idx[10];
  unsigned int prob_max_size = 0;
  prob_max_idx[0] = idx;
  prob_max_size = 1;
  for (;*cp != '\n'; ++cp)
    {
      idx = ((idx << 2) & kmask) + g_idx_dict[*cp];
      prob = calc_prob(idx, k, prf);

      if (prob_max < prob)
	{
	  prob_max = prob;
	  prob_max_idx[0] = idx;
	  prob_max_size = 1;
	}
      else if (prob_max == prob)
	{
	  prob_max_idx[prob_max_size] = idx;
	  ++prob_max_size;
	}
    }

  for (i = 0; i < prob_max_size; ++i)
    print_clump(prob_max_idx[i], k);
  return 0;
}
