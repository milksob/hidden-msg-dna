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

const char * find_line_start(char *buf_start, const char * linep)
{
  const char *p = linep;
  while (*p != '\n')
    {
      --p;
      if (buf_start == p)
	{
	  printf("Wrong file format\n");
	  exit(0);
	}
    }
  ++p;
  return p;
}

char g_idx_dict['T' + 1];
typedef enum Base { A, C, G, T, KMER_SIZE = 2 } Base;
void init_dict()
{
  g_idx_dict['A'] = A;
  g_idx_dict['C'] = C;
  g_idx_dict['G'] = G;
  g_idx_dict['T'] = T;
}

int calc_idx(int idx, char ch)
{
  return (idx << 2) + g_idx_dict[ch];
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

int calc_distance(int kmer_idx, int pattern_idx, int k)
{
  static int diff;
  diff = kmer_idx ^ pattern_idx;
  static const int odd_mask = 0xAAAAAAAA;   //10 10 ... mask
  static const int even_mask = 0x55555555;  //01 01 ... mask

  diff |= ((diff & odd_mask) >> 1); // OR evens with odds
  diff &= even_mask; // clear all odds;
  diff = __builtin_popcount(diff);  // number of different letters

  return diff;
}

enum { POS_SIZE = 5, BUF_SIZE = 5*1024*1024 };

int main(int argc, char **argv)
{
  int fd = open_file(argc, argv);

  int count = 0;
  //  char *output = (char *)malloc(BUF_SIZE);
  char *buf = (char *)malloc(BUF_SIZE);
  int len = read_file(fd, buf, BUF_SIZE);
  const char *buf_end = buf + len;

  int d = 0;
  const char *cp = find_line_start(buf, buf_end - 2);
  sscanf(cp, "%d", &d);

  init_dict();
  int k = 0;
  cp = buf;
  int pattern_idx = 0;
  while (*cp != '\n')
    {
      pattern_idx = calc_idx(pattern_idx, *cp);
      ++cp;
    }
  k = cp - buf;
  assert(k*KMER_SIZE < sizeof(pattern_idx)*8);

  ++cp; // pass eol

  // fill first kmer
  int k_mask = (1 << (2*k)) - 1;
  int kmer_idx = 0;
  int i;
  for (i = 0; i < k; ++i)
    {
      kmer_idx = calc_idx(kmer_idx, *cp);
      ++cp;
    }
  int d_idx = k;
  d_idx = calc_distance(kmer_idx, pattern_idx, k);
  int pos = 0; // first kmer pos is 0
  if (d_idx <= d)
    ++count;
  ++pos;

  while (*cp != '\n')
    {
      kmer_idx = ((kmer_idx << 2) & k_mask) + g_idx_dict[*cp];
      d_idx = calc_distance(kmer_idx, pattern_idx, k);
      if (d_idx <= d)
	++count;
      ++cp;
      ++pos;
    }

  printf("Output : %d\n", count);

  free(buf);
  return 0;
}
