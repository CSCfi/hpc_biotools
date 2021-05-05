#!/usr/bin/perl -w
use strict;

my $arg=\@ARGV;

my $t1=$arg->[0];
my $t2=$arg->[1];
my $t3=$arg->[2];

open(INPUT,"$t1");
open(INPUT2,"$t2");
while (my $rivi=<INPUT>){
 my $rivi2=<INPUT2>;
 chomp $rivi;
 chomp $rivi2;
 if ( $t3 eq "TAB" ){
   print "$rivi\t\"$rivi2\"\n";
 } else {
   print "$rivi\,\"$rivi2\"\n";
 }
}

