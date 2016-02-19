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

enum { POS_SIZE = 5, BUF_SIZE = 5*1024*1024 };

int main(int argc, char **argv)
{
  int fd = open_file(argc, argv);

  char buf[BUF_SIZE];
  int len = read_file(fd, buf, BUF_SIZE);

  int hamming = 0;

  char *cpl = buf;
  char *cpr = buf + len/2 + len%2;
  cpr = find_line_start(cpl, cpr);

  while (*cpr && (*cpl) != '\n' && (*cpr) != '\n')
    {
      if (*cpr != *cpl)
	{
	  ++hamming;
	}
      ++cpl;
      ++cpr;
    }

  printf("Output: hamming distance is %d\n", hamming);
  return 0;
}
