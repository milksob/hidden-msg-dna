#include <stdio.h>
#include <fcntl.h>

enum { BUF_SIZE = 1024*1024 };
int main(int argc, char *argv[])
{
  if (argc < 2)
    {
      printf("Usage: %s filename\n", argv[0]);
      return 0;
    }
  
  int fd = open(argv[1], O_RDONLY);
  if (fd == -1)
    printf("Can't open file %s\n", argv[1]);

  char buf[BUF_SIZE];
  int length = read(fd, buf, BUF_SIZE);
  if (length == -1)
    {
      printf("Error reading file\n");
      return;
    }
  if (length > BUF_SIZE)
    {
      printf("File is too big\n");
      return;
    }

  char *pos = buf;
  char *pos_end = pos + length;

  printf("text: ");
  while (pos < pos_end)
    {
      ++pos;
      if (*pos == '\n')
	break;
      printf("%c", *pos);
    }
  printf("\n");

 // next symbol. format text \n pattern
  if (*pos != '\n')
    {
      printf("wrong file format\n");
      return;
    }
  char *pat = ++pos;
  int pat_l = pos_end - pat;
  pos = pos_end - 1; // last symbol
  if (*pos == '\n')
    {
      --pat_l;
      *pos = 0;
    }
  printf("pattern = %s, length = %d\n", pat, pat_l);

  int count = 0;
  pos = buf;
  pos_end = pat - pat_l;
  while (pos < pos_end)
    {
      if (strncmp(pos, pat, pat_l) == 0)
	++count;
      ++pos;
    }
  printf("Output: %d\n", count);
}
