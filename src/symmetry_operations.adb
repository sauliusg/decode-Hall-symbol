with Ada.Text_IO; use Ada.Text_IO;

package body Symmetry_Operations is
   
   procedure Put (F : File_Type; S : Symmetry_Operator) is
   begin
      for J in Symmetry_Operator'Range(1) loop
         for K in Symmetry_Operator'Range(1) loop
            Put (F, " " & Float'Image (S (J,K)));
         end loop;
         New_Line (F);
      end loop;
   end;
   
   procedure Put (F : File_Type; S : Crystallographic_Translation) is
      Ratio : Float;
   begin
      for I in S'Range loop
         Ratio := Float (S (I).Numerator) / Float (S (I).Denominator); 
         Put (F, " " & Float'Image (Ratio));
      end loop;
      New_Line (F);
   end;
   
   function To_Symmetry_Operator (T : Crystallographic_Translation) 
                                 return Symmetry_Operator
   is
      S : Symmetry_Operator := Unity_Matrix;
   begin
      for I in T'Range loop
         S (I,4) := Float (T (I).Numerator) / Float (T (I).Denominator);
      end loop;
      return S;
   end;
   
   function Axis_Index (Direction : Known_Axis_Direction) return Positive is
   begin
      case Direction is
         when X_AXIS => return 1;
         when Y_AXIS => return 2;
         when Z_AXIS => return 3;
      end case;
   end Axis_Index;
      
   function To_Symmetry_Operator (T : Crystallographic_Translation_Component;
                                  Axis_Direction : Known_Axis_Direction)
                                 return Symmetry_Operator
   is
      S : Symmetry_Operator := Unity_Matrix;
   begin
      S (Axis_Index (Axis_Direction), 4) :=
        Float (T.Numerator) / Float (T.Denominator);
      return S;
   end To_Symmetry_Operator;
   
   procedure Skip_Spaces (S : in String; Pos : in out Integer ) is
   begin
      while Pos <= S'Last and then S (Pos) = ' ' loop
         Pos := Pos + 1;
      end loop;
   end;
   
   procedure Decode_Centering_Symbol
     (
      Symbol : in String;
      Pos : in out Positive;
      Centering : out Symmetry_Operator_Array;
      N_Centering : out Positive
     )
   is
   begin
      Skip_Spaces (Symbol, Pos);
      Centering (1) := Unity_Matrix;
      case Symbol (Pos) is
         when 'P' =>
           N_Centering := 1;
         when 'A' =>
           Centering (2) := To_Symmetry_Operator (A_Translation_Vector);
           N_Centering := 2;
         when 'B' =>
           Centering (2) := To_Symmetry_Operator (B_Translation_Vector);
           N_Centering := 2;
         when 'C' =>
           Centering (2) := To_Symmetry_Operator (C_Translation_Vector);
           N_Centering := 2;
         when 'I' =>
           Centering (2) := To_Symmetry_Operator (I_Translation_Vector);
           N_Centering := 2;
         when 'F' =>
           Centering (2) := To_Symmetry_Operator (F_Translation_Vector_1);
           Centering (3) := To_Symmetry_Operator (F_Translation_Vector_2);
           Centering (4) := To_Symmetry_Operator (F_Translation_Vector_3);
           N_Centering := 4;
         when 'R' =>
           Centering (2) := To_Symmetry_Operator (R_Translation_Vector_1);
           Centering (3) := To_Symmetry_Operator (R_Translation_Vector_2);
           N_Centering := 3;
         when others =>
            raise UNKNOWN_CENTERING with
              "unknown centering symbol " & Character'Image(Symbol (Pos));
      end case;
      Pos := Pos + 1;
   end;
   
end Symmetry_Operations;
