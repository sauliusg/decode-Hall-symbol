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
   
end Symmetry_Operations;
