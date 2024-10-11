gcc main.c -Wall -lm -o out
.\out.exe
magick mogrify -format jpg -background white -alpha remove -quality 97 out/p*.eps
ffmpeg -framerate 15 -y -i .\out\p%%05d.jpg -c:v libx264 .\out\_movie.avi