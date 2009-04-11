#!/usr/bin/perl -w -i.bak -p
if( ! /^\#/) {
    s/Prom\s+(\-\d+|0|\.\.\.)/Prom    1/;
}
