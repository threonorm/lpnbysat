let genere_poly s =
  for i=0 to (String.length s)-1 do
    if (int_of_char(s.[i])-48==1) 
    then
      Format.printf "+ a%d " i 
  done;
  Format.printf "\n" 

let generate_var s =
  for i= 0 to (String.length s)-1 do
    Format.printf ",a%d " i;
  done;
  Format.printf "=var('";
  for i= 0 to (String.length s)-1 do
    Format.printf "a%d " i;
  done;
  Format.printf "')\n"

let generate_square n =
  for i=1 to n/2 do
  Format.printf "e%d=eq%d*eq%d;\n" i (2*i-1) (2*i)
  done;
  Format.printf "[";
  for i=1 to n/2 - 1 do
  Format.printf "e%d, " i;
  done;
  Format.printf "e%d]\n" (n/2) 
       
  

let ()=
      generate_square (int_of_string (Sys.argv).(2));
      generate_var ((Sys.argv).(1)); 
      genere_poly ((Sys.argv).(1)) 
