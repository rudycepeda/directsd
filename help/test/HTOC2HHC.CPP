#include<stdio.h>
#include<string.h>

FILE *in, *out;

void puttab ( int n, char s[] )
{
int i;
for ( i=1; i <= n; i++ ) fputc( '\t', out );
fprintf(out, "%s", s);
fputc('\n', out);
}

void main(int argc, char *argv[])
{
int i, nest = 0;
char s[1000], w[1000], title[200], fname[100];
char *p;

if ( argc < 2 ) {
   puts("Usage: htoc2hhc <filename>\n");
   return;
   }
in = fopen ( argv[1], "r" );
if ( in == NULL ) {
   printf("File %s not found\n", argv[1]);
   return;
   }

out = fopen ( "hhc", "w" );
while ( fgets(s, 1000, in) )
   {
   i = 0;
   while ( s[i] == '>' ) i++;
   if ( i > nest ) puttab ( nest, "<UL>");
   if ( i < nest ) puttab ( i
   , "</UL>");
   p = strchr ( s, '\n' );
   *p = 0;
   p = strchr ( s, '&' );
   strcpy ( title, p+1 );
   *p = 0;
   p = strstr ( s, ".htt" );
   strcpy ( p, ".html" );
   puttab ( i,   "<LI> <OBJECT type=\"text/sitemap\">" );
   sprintf(w, "<param name=\"Name\" value=\"%s\">", title);
   puttab ( i+1, w );
   sprintf(w, "<param name=\"Local\" value=\"%s\">", s );
   puttab ( i+1, w );
   puttab ( i, "</OBJECT>");
   nest = i;
   }

fclose(in);
fclose(out);
}