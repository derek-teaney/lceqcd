#!/usr/bin/perl -w 

# Open a fortran file look for includes and use statements
#
# USAGE:
#
# fdepends.pl name.f name.d "-I dir1/ -I dir2/"
$fname  = $ARGV[0] ;

$gIncludePath = $ARGV[2];
$gOutputOpt = $ARGV[3] ;
@dirs = split(/\-I/,$gIncludePath) ;

#print "@dirs\n";

if ( $fname =~m/(.*)\/(.*)/ ){
   $name = $2 ; 
   $path = "$1/" ;
}else{
   $name = $fname ;
   $path = "" ;
}

$outname = $ARGV[1] ;
($dependo=$outname)=~ s/\.d\s*$/\.o/ ;

open(FI,"<$ARGV[0]") ;
open(FO,">$ARGV[1]") ;

%MODULESIN=() ;
%INCLUDES=() ;
%USEMODUL=() ;
while (<FI>) {
   $line = $_ ;
	if($line=~m/^\s*module\s+([^!]+)/i){
      #print $line ;
       ($tmp = lc("$1")) =~s/\s//g ;
       $modname = &fixmodname($tmp) ;
       $MODULESIN{$modname} = $modname;
   }
	if ($line=~m/^\s*include\s+["\']([^"\']+)["\']/i){
       $tmp = "$1";
       $INCLUDES{$tmp} = $tmp;
	} 
   if ($line=~/^\s*use\s+([^,!]+)/i){
       #print $line ;
       ($tmp = lc("$1")) =~s/\s//g;
       $modname = &fixmodname($tmp) ;
       $USEMODUL{$modname} = $modname;
	}
}
foreach $key (keys %MODULESIN) {
   if (exists($USEMODUL{$key})) {
      delete($USEMODUL{$key}) ; 
   }
}

$dependline = "$dependo" ;
$dependline .= ": $fname" ;
foreach $key (keys %USEMODUL) {
   $pathkey = &findinpath2($key) ;
   $dependline .= " $pathkey"  ;
}
foreach $key (keys %INCLUDES) {
   $pathkey = &findinpath($key) ;
   $dependline .= " $pathkey"  ;
}
print FO "$dependline \n" ;
#print FO "\t$cmd\n" ;

foreach $key (keys %MODULESIN) {
   $dependline = "$key"  ;
   $dependline .= ": $dependo "; 
   print FO "$dependline \n" ;
}


# Converts the module name to lower case
# and apends .mod

sub fixmodname
{
   my($mod) = @_ ;
   # converts name to lower case
   if ($gOutputOpt eq "ifc") {
      $mod = uc($mod) ;
      $mod = $mod.".mod"  ;
   }else{
      $mod = lc($mod) ;
      $mod = $mod.".mod"  ;
   }
   # appends .mod
   return $mod ;
}


# 
# 
sub findinpath
{
   my($key) = @_ ;
   my($dir, $pathname) ;
   #print "key=$key.\n" ;
   DIR: foreach $dir(@dirs) {
      #print  "$dir\n" ;
      $dir =~s/\s*//g ;
      if ($dir eq ""){ next DIR ; }

      $pathname = $dir.$key ;
      #print "pathname=$pathname \n" ;
      if (-e $pathname) {
         return  $pathname ;
      } 
   }
   print "Cant find key = $key\n" ;
   return "" ;
}

# 
# 
sub findinpath2
{
   my($key) = @_ ;
   my($dir, $pathname) ;
   #print "key=$key.\n" ;
   DIR: foreach $dir(@dirs) {
      #print  "$dir\n" ;
      $dir =~s/\s*//g ;
      if ($dir eq ""){ next DIR ; }

      $pathname = $dir.$key ;
      #print "pathname=$pathname \n" ;
      if (-e $pathname) {
         return  $pathname ;
      } 
   }
   return $key
}


