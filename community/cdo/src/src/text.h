#ifndef _TEXT_H
#define _TEXT_H

#define RESET		0
#define BRIGHT 		1
#define DIM		2
#define UNDERLINE 	4
#define BLINK		5
#define REVERSE		7
#define HIDDEN		8

#define BLACK 		0
#define RED		1
#define GREEN		2
#define YELLOW		3
#define BLUE		4
#define MAGENTA		5
#define CYAN		6
#define	WHITE		7


#define COLOR_STDOUT (stdout_is_tty && CDO_Color)
#define COLOR_STDERR (stderr_is_tty && CDO_Color)

void set_text_color(FILE *fp, int attr, int fg);
void reset_text_color(FILE *fp);

#endif  /* _TEXT_H */
