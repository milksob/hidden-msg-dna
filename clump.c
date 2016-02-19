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

typedef struct Lq_node Lq_node;

typedef struct Lq_node
{
  Lq_node *next;
  
  int kmer;
} Lq_node;

int main(int argc, char **argv)
{
  int fd = open_file(argc, argv);

  char buf[BUF_SIZE];
  int len = read_file(fd, buf, sizeof(buf));

  char *vars = find_line_start(buf, buf + len - 2);
  int k, L, t;
  sscanf(vars, "%d %d %d", &k, &L, &t);
  --vars;
  *vars = '\0';
  //  char *genome = buf;
  size_t tmap_size = pow(4, k) * sizeof(int);

  Lq_node *nodes = (Lq_node *)malloc(sizeof(Lq_node)*L);
  int *tmap = (int *)malloc(tmap_size);
  memset(tmap, 0, tmap_size);

  printf("k = %d, L = %d, t = %d\n", k, L, t);

  int idx = 0;
  char *cp = buf;
  int i;
  for (i = 0; i < k; ++i)
    {
      idx = (idx << 2) + calc(*cp);
      // most recent nucleo is ls-2bits
      ++cp;
    }
  ++tmap[idx]; // first idx done

  Lq_node *front = nodes;
  front->kmer = idx;

  Lq_node *np;
  int k_mask = (1 << (2*k)) - 1;

  int nj = 1;
  for (; i < L; ++i, ++nj, ++cp)
    {
      np = &nodes[nj];
      nodes[nj-1].next = np;

      idx = ((idx << 2) & k_mask) + calc(*cp);
      np->kmer = idx;

      if (tmap[idx] < t)
	{
	  ++tmap[idx];
	  if (tmap[idx] >= t)
	    print_clump(idx, k);
	}
    }

  Lq_node *last = np;
  while (*cp)
    {
      np = front;
      front = front->next;

      if (tmap[np->kmer] < t)
	--tmap[np->kmer];

      idx = ((idx << 2) & k_mask) + calc(*cp);
      np->kmer = idx;

      if (tmap[idx] < t)
	{
	  ++tmap[idx];
	  if (tmap[idx] >= t)
	    print_clump(idx, k);
	}

      last->next = np;
      last = np;

      ++cp;
    }

  free(tmap);
  free(nodes);
  return 0;
}
