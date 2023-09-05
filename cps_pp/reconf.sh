#!/usr/bin/bash
autoreconf -fvi
# configuring version number
echo "`git log -n 1 --format=format:"#define GITHASH \\"%H:%d\\"%n" HEAD`" > include/version.h
