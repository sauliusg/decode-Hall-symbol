pragma Ada_2022;

package GCD_Mod with Spark_Mode Is

   function Equivalent (L1, L2 : Boolean) return Boolean
   is ((L1 and then L2) or else (not L1 and then not L2))
     with
     Spark_Mode,
     Ghost,
     Pre => L1 in Boolean and L2 in Boolean,
     Post => Equivalent'Result = ((L1 and L2) or (not L1 and not L2));
   
   function Is_Divisor (A : Natural; D : Positive) return Boolean
     is (A mod D = 0)
     with
       Ghost,
       Pre => A in Natural and then D in Positive,
       Post => Is_Divisor'Result = (A mod D = 0);
   
   function Is_GCD (A, B, D : in Natural) return Boolean
   is 
      (A mod D = 0 and then
       B mod D = 0 and then
         (for all D1 in D .. Positive'Last =>
            D = D1 or else
            (A mod D1 /= 0) or else
            (B mod D1 /= 0)))
        with
        Ghost,
        Pre => A >= 0 and then B > 0 and then D > 0,
        Post =>
        Is_GCD'Result = 
          (A mod D = 0 and then
           B mod D = 0 and then
             (for all D1 in D .. Positive'Last =>
                D = D1 or else
                (A mod D1 /= 0) or else
                (B mod D1 /= 0)))
          ;
   
   function GCD (A, B : in Positive) return Positive
     with
     Pre  => 
       A > 0 and then A <= Positive'Last and then 
       B > 0 and then B <= Positive'Last,
     Post => Is_GCD (A, B, GCD'Result);
         
end GCD_Mod;
