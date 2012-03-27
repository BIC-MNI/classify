#! /bin/sh

set -e

# add libtool stuff
libtoolize --automake

# normal automake/conf bits
aclocal -I m4
autoheader
automake --add-missing --copy
autoconf

