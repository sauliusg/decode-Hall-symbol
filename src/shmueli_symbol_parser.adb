package body Shmueli_Symbol_Parser is
   
   function Decode_Shmueli_Symbol (Symbol : in String)
                                   return Symmetry_Operator_Array
   is
      Max_Symmetry_Operators : constant Integer := 192;
      
      Symmetry_Operators :
        Symmetry_Operator_Array (1 .. Max_Symmetry_Operators);
      N_Symmetry_Operators : Positive := 1;
   begin
      Symmetry_Operators (1) := Unity_Matrix;
      
      return Symmetry_Operators (1..N_Symmetry_Operators);
   end;
   
end Shmueli_Symbol_Parser;
