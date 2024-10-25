pragma Ada_2022;

package body GCD_Mod with Spark_Mode Is
   
   function GCD (A, B : in Positive) return Positive
   is
      X, Y, T : Natural;
   begin
      X := A;
      Y := B;
      
      pragma Assert
        (for all N in Positive => Is_GCD (0, N, N));
           
      pragma Assume
        (for all M in Positive =>
           (for all N in Positive =>
              (for all D in Positive =>
                 Equivalent (Is_GCD (M, N, D), Is_GCD (N mod M, M, D)) )));
      
      while X > 0 loop
         
         pragma Loop_Invariant
           (for all G in Positive =>
              Equivalent (Is_GCD(A, B, G), Is_GCD(X, Y, G)));
         
         pragma Loop_Invariant (Y > 0);
              
         T := X;
         X := Y mod X;
         Y := T;
            
      end loop;
      
      pragma Assert (Y > 0);
      pragma Assert (Is_GCD (A, B, Y));
      
      return Y;
   end;
   
end GCD_Mod;
