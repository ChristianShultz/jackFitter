#!/bin/tcsh  

if( $#argv != 1 ) then
  echo "usage <1=headders, 2=source>"
  exit 1
endif  

set headders = `ls jackFitter/*.h`
set source = `ls jackFitter/*.cc`


if( $1 == "1" ) then 
  set files = "$headders"
endif

if( $1 == "2" ) then 
  set files = "$source"
endif

foreach x ( $files ) 
  echo "$x    \\"
end
