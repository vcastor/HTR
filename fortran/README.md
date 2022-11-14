# Welcome to Hartree-Fock-Roothaan (FORTRAN version) :rocket:

For more knowlagment about HF method, there is a pdf with some explanations
about the theory level. Don't be shy to send me a message if something isn't
written propertly or you have questions.

## Compiling

Just type: `make` on your terminal at the directory where is the source code.

If that didn't work correctly I'll give you the answer of the problems that you
might have.

First of all you need a fortran compiler, I recomend gfortran, even if its
written on C/C++, is more robust and extended use than other compilers. Intel
Fortran Compiler, NVIDIA nvfortran and AMD flang are developed just
specifically for their owns architectures.

Also gfortran is free and can be used on Linux, macOS and its which we used
:grimacing:. However, of course you can change it for your favorite compiler in
the `Makefile`. If you don't have any compiler, once again, I recomend [GNU
Compiler Collection] (http://...) (its not only fortran compiler, it can also
compile: `C`, `C++`, `Objective-C`, `Ada`, `Go`, and `D`).

<hr style="border:2px solid gray"> </hr>

The program was developed using
[LAPACK](http://www.netlib.org/lapack/#_lapack_version_3_10_0) library, if your
computer doesn't have the library you will have an error like:

>`"_dsygvd_", referenced from:`
> `     _MAIN__ in ccZlxZJJ.o`

Therefore, the problem above can be fixed just intalling LAPACK library. In
case that you have it but doesn't compailing neither, is because you have it in
a unusual directory, not in one that I'm supposing in the Makefile, just modify
the `Makefile` at the line where we call the library in the compaling.

<hr style="border:2px solid gray"> </hr>

Fundamental Theorem of Computing: turn off, turn on and try again.
If you alredy installed LAPACK and you're trying to compile on the same terminal
session maybe its necessary to restart the `$PATH`, you can type the
command `~/.bash_profile` or just open a new terminal.

If you have any other issue, you can send me a message, I'll try to help you, but
also try to google it, there aren't many ways that the compilation doesn't work
propertly.

## How run a calculation :wink:

On the same directory is a laucher call it `launcher.sh` (I have a lot of
creativity, I know :blush:). The laucnher would have the excecutable permissions, if it
doesn't have give it with the next command:

> `chmod +x launcher.sh`

To send the calculation just give the input as an argument of the launcher.
Example:

> `$ ./launcher.sh H2.inp &`

The `&` just to not freez the terminal, its doesn't matter the name, but the
extension should be `inp`, that becuase the launcher will create the out put
file with the same name but with the extention `out` and a scratch file as
`moreout`.

## To see the source code :beers:

If you want to see the source code, take acount that will be easer to read if
you have a code editor as Sublime Text, Atom, Visual Studio Code, ... or just
vi/vim with **syn on**, that to have highlighted syntaxes on reserved words.

