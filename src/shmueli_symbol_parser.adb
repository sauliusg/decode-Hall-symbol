with Ada.Text_IO; use Ada.Text_IO;
with Shmueli_Matrices; use Shmueli_Matrices;

package body Shmueli_Symbol_Parser is
   
   function As_Symetry_Operator (R : Rotation_Matrix) 
                                return Symmetry_Operator
   is
      S : Symmetry_Operator;
   begin
      for I in 1 .. 3 loop
         for J in 1 .. 3 loop
            S (I, J) := Float (R (I, J));
         end loop;         
      end loop;
      for I in 1 .. 3 loop
         S (I, 4) := 0.0;
         S (4, I) := 0.0;
      end loop;
      S (4, 4) := 1.0;
      return S;
   end;
   
   procedure Get_Shmueli_Symbol_Rotation
     (
      Symbol : in String;
      Pos : in out Positive;
      Rotation : out Symmetry_Operator
     )
   is
   begin
      if Symbol (Pos .. Pos + 1) = "1A" then
         Rotation := As_Symetry_Operator (M_1A);
      elsif Symbol (Pos .. Pos + 1) = "2A" then
         Rotation := As_Symetry_Operator (M_2A);
      elsif Symbol (Pos .. Pos + 1) = "2B" then
         Rotation := As_Symetry_Operator (M_2B);
      elsif Symbol (Pos .. Pos + 1) = "2C" then
         Rotation := As_Symetry_Operator (M_2C);
      elsif Symbol (Pos .. Pos + 1) = "2D" then
         Rotation := As_Symetry_Operator (M_2D);
      elsif Symbol (Pos .. Pos + 1) = "2E" then
         Rotation := As_Symetry_Operator (M_2E);
      elsif Symbol (Pos .. Pos + 1) = "2F" then
         Rotation := As_Symetry_Operator (M_2F);
      elsif Symbol (Pos .. Pos + 1) = "2G" then
         Rotation := As_Symetry_Operator (M_2G);
      elsif Symbol (Pos .. Pos + 1) = "3Q" then
         Rotation := As_Symetry_Operator (M_3Q);
      elsif Symbol (Pos .. Pos + 1) = "3C" then
         Rotation := As_Symetry_Operator (M_3C);
      elsif Symbol (Pos .. Pos + 1) = "4C" then
         Rotation := As_Symetry_Operator (M_4C);
      elsif Symbol (Pos .. Pos + 1) = "6C" then
         Rotation := As_Symetry_Operator (M_6C);
      else
         raise INVALID_SYMBOL with
           "rotation symbol pair '" & Symbol (Pos .. Pos + 1) & "' " &
           "unexpected at position " & Pos'Image & " in '" &
           Symbol & "'";
      end if;
      Pos := Pos + 2;
   end;
   
   function Decode_Shmueli_Symbol (Symbol : in String)
                                   return Symmetry_Operator_Array
   is
      Pos : Positive := 1; -- current position in the string 'Symbol'.
      
      Max_Symmetry_Operators : constant Integer := 192;
      
      Symmetry_Operators :
        Symmetry_Operator_Array (1 .. Max_Symmetry_Operators);
      N_Symmetry_Operators : Positive := 1;
      
      Crystal_System : Character;
      
      Inversions : array (1..2) of Symmetry_Operator :=
        (Unity_Matrix, Ci_Matrix);
      N_Inversions : Positive;
      
      Max_Centering : constant Positive := 8;
      Centering : Symmetry_Operator_Array (1..Max_Centering);
      N_Centering : Positive;
      
   begin
      Symmetry_Operators (1) := Unity_Matrix;
      
      Decode_Centering_Symbol (Symbol, Pos, Centering, N_Centering);
      
      Crystal_System := Symbol (Pos);
      Pos := Pos + 1;
      
      case Symbol (Pos) is
         when 'N' => N_Inversions := 1;
         when 'C' => N_Inversions := 2;
         when others =>
            raise INVALID_SYMBOL with
              "centrosymmetry symbol '" & Symbol (Pos .. Pos) & "' " &
              "is not expected, expecting 'N' or 'C' at position " &
              Pos'Image & " in '" & Symbol & "'";
      end case;
      Pos := Pos + 1;
      
      for I in 1 .. 3 loop
         exit when Pos > Symbol'Length;
         
         if Symbol (Pos) /= '$' then
            raise INVALID_SYMBOL with
              "expecting '$' instead of " &
              "instead of '" & Symbol (Pos .. Pos) & "' " &
              "at position " & Pos'Image & " in '" &
              Symbol & "'";
         end if;
         
         Pos := Pos + 1;
         
         declare 
            Rotation_Type : Character := Symbol (Pos);
            Rotation : Symmetry_Operator;
            Translation : Crystallographic_Translation;
         begin
            Pos := Pos + 1;
            Get_Shmueli_Symbol_Rotation (Symbol, Pos, Rotation);
            case Rotation_Type is
               when 'P' => null;
               when 'I' => Rotation := Ci_Matrix * Rotation;
               when others =>
            raise INVALID_SYMBOL with
              "rotation matrix symbol '" & Symbol (Pos-3 .. Pos-3) & "' " &
              "is not expected, expecting 'P' (Proper) or 'I' (Improper) " &
              "at position " &
              Pos'Image & " in '" & Symbol & "'";
            end case;
            Decode_Shmueli_Symbol_Translation
              (
               Symbol (Pos .. Pos + 2),
               Translation
              );
            Pos := Pos + 3;
            Add (Rotation, Translation);
            if Rotation /= Unity_Matrix then
               N_Symmetry_Operators := N_Symmetry_Operators + 1;
               Symmetry_Operators (N_Symmetry_Operators) := Rotation;
            end if;
         end;
      end loop;
      
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
      
      -- Reconstruct all symmetry operators:
      
      Build_Group (Symmetry_Operators, N_Symmetry_Operators);
      
      return Symmetry_Operators (1 .. N_Symmetry_Operators);
   end;
   
end Shmueli_Symbol_Parser;
