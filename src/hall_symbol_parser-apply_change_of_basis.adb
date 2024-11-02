separate (Hall_Symbol_Parser)

procedure Apply_Change_Of_Basis
  (
   Symmetry_Operators : in out Symmetry_Operator_Array;
   N_Symmetry_Operators : in out Natural;

   Centering : in out Symmetry_Operator_Array;
   N_Centering : in out Natural;

   Inversions : in out Symmetry_Operator_Array;
   N_Inversions : in out Natural;

   Change_Of_Basis : in Symmetry_Operator
  )
is
begin
   
   if Change_Of_Basis /= Unity_Matrix then
      declare
         V : Symmetry_Operator := Change_Of_Basis;
         V_Inv : Symmetry_Operator;
      begin
         V_Inv := Invert (V);
         
         for I in 2..N_Symmetry_Operators loop
            Symmetry_Operators (I) := V * Symmetry_Operators (I) * V_Inv;
         end loop;
         
         for I in 2..N_Centering loop
            Centering (I) := V * Centering (I) * V_Inv;
         end loop;
         
         while N_Centering > 1 and then 
           Centering (N_Centering) = Centering (1) 
         loop
            -- remove centerings that became unit matrices after C-o-B:
            N_Centering := N_Centering - 1;
         end loop;
         
         if N_Inversions = 2 then
            Inversions (2) := V * Inversions (2) * V_Inv;
         end if;
      end;
   end if;
   
   -- Generate additional centering operators:
   
   declare
      type Vector_Type is array (1..4) of Float;
      
      function "*" (S : Symmetry_Operator; T : Vector_Type)
                   return Vector_Type 
      is
         R : Vector_Type := (others => 0.0);
      begin
         for I in R'Range loop
            for K in T'Range loop
               R (I) := R (I) +
                 S (I,K) * T (K);
            end loop;
         end loop;
         return R;
      end;
      
      function To_Symmetry_Operator (T : Vector_Type)
                                    return Symmetry_Operator
      is
         S : Symmetry_Operator := Unity_Matrix;
      begin
         for I in 1..3 loop
            S (I,4) := T (I);
         end loop;
         Snap_To_Crystallographic_Translations (S);
         return S;
      end;
      
      function Is_Centering (T : Vector_Type) return Boolean is            
         function Fract (X : Float) return Float is (X - Float'Floor (X));
         begin
            for Component of T loop
               if abs (Fract (Component)) >= Eps then
                  return True;
               end if;
            end loop;
            return False;
         end;
         
         Unit_Vectors : array (1..3) of Vector_Type :=
           (
            (1.0, 0.0, 0.0, 1.0),
            (0.0, 1.0, 0.0, 1.0),
            (0.0, 0.0, 1.0, 1.0)
           );
         
         C_O_B_Rotation : Symmetry_Operator := Change_Of_Basis;
         
      begin
         for I in 1..3 loop
            C_O_B_Rotation (I,4) := 0.0;
         end loop;
         for Vector of Unit_Vectors loop
            if Is_Centering (C_O_B_Rotation * Vector) then
               N_Centering := N_Centering + 1;
               Centering (N_Centering) :=
                 To_Symmetry_Operator (C_O_B_Rotation * Vector);
            end if;
         end loop;
         
         -- see if multiplication of two new centerings gives a third one:
         Build_Group (Centering, N_Centering);
      end;
      
      -- Reconstruct all rotation operators:
      
      Build_Group (Symmetry_Operators, N_Symmetry_Operators);
      
      -- Add centering and inversion matrices:
      
      declare
         M : Positive := N_Symmetry_Operators;
         New_Symmetry_Operator : Symmetry_Operator;
      begin
         for I in 1..N_Inversions loop
            for C in 1..N_Centering loop
               if I /= 1 or else C /= 1 then
                  for S in 1..N_Symmetry_Operators loop
                     New_Symmetry_Operator :=
                       Symmetry_Operators (S) * Centering (C) * Inversions (I);
                     if not Has_Symmetry_Operator (Symmetry_Operators, M, 
                                                   New_Symmetry_Operator) then
                        M := M + 1;
                        Symmetry_Operators (M) := New_Symmetry_Operator;
                     end if;
                  end loop;
               end if;
            end loop;
         end loop;
         N_Symmetry_Operators := M;
      end;      
      
end Apply_Change_Of_Basis;
