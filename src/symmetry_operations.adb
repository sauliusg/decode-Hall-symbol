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
   
   function As_String (S : Symmetry_Operator) return String is
      Buffer : String (1..100); -- large enough to hold any symop
      Pos : Positive := 1;
      Non_Zero_Printed : Boolean;
      
      function Rational_Translation (T : Float) return String is
         Buffer : String := Float'Image (T);
         Idx : Positive := 1;
      begin
         if T = 0.0 then
            return "";
         elsif T = 0.5 then
            return "+1/2";
         elsif T = 1.5 then
            return "+3/2";
         elsif T = 1.0/3.0 then
            return "+1/3";
         elsif T = 2.0/3.0 then
            return "+2/3";
         elsif T = 1.0/4.0 then
            return "+1/4";
         elsif T = 3.0/4.0 then
            return "+3/4";
         elsif T = 1.0/6.0 then
            return "+1/6";
         elsif T = 5.0/6.0 then
            return "+5/6";
         elsif T = 1.0/8.0 then
            return "+1/8";
         elsif T = 3.0/8.0 then
            return "+3/8";
         elsif T = 5.0/8.0 then
            return "+5/8";
         elsif T = 7.0/8.0 then
            return "+7/8";
         else
            while Idx <= Buffer'Last and then
              (Buffer (Idx) = ' ' or else
                 Buffer (Idx) = '+' or else
                 Buffer (Idx) = '-') loop
               Idx := Idx + 1;
            end loop;
            return (if T < 0.0 then "-" else "+") & Buffer (Idx..Buffer'Last);
         end if;
      end Rational_Translation;
      
   begin
      for I in 1 .. S'Last(2) - 1 loop
         Non_Zero_Printed := False;
         for J in S'Range(1) loop
            if J < 4 then
               -- rotation part:
               if S (I,J) /= 0.0 then
                  if S (I,J) < 0.0 then
                     Buffer (Pos) := '-';
                     Pos := Pos + 1;
                  else
                     if Non_Zero_Printed then
                        Buffer (Pos) := '+';
                        Pos := Pos + 1;
                     end if;
                  end if;
                  if Abs (S (I,J)) /= 1.0 then
                     declare
                        R : String := Rational_Translation (abs(S (I,J)));
                     begin
                        for Ch of R (2..R'Last) loop
                           Buffer (Pos) := Ch;
                           Pos := Pos + 1;
                        end loop;
                     end;
                     Buffer (Pos) := '*';
                     Pos := Pos + 1;
                  end if;
                  case J is
                     when 1 => Buffer (Pos) := 'X';
                     when 2 => Buffer (Pos) := 'Y';
                     when 3 => Buffer (Pos) := 'Z';
                     when others =>
                        raise CONSTRAINT_ERROR;
                  end case;
                  Pos := Pos + 1;
                  Non_Zero_Printed := True;
               end if;
            else
               -- translation part:
               for C of Rational_Translation (S (I,J)) loop
                  Buffer (Pos) := C;
                  Pos := Pos + 1;
               end loop;
            end if;
         end loop;
         if I < 3 then
            Buffer (Pos) := ',';
            Pos := Pos + 1;
         end if;
      end loop;
      
      return Buffer (1..Pos-1);
   end As_String;
   
end Symmetry_Operations;
