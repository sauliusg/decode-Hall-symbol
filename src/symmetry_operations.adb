with Ada.Text_IO; use Ada.Text_IO;

package body Symmetry_Operations is
   
   function Machine_Epsilon return Float is
      Epsilon : Float := 1.0;
   begin
      while 1.0 + Epsilon /= 1.0 loop
         Epsilon := Epsilon / 2.0;
      end loop;
      return Epsilon;
   end;
   
   function Eps return Float is (16.0 * Machine_Epsilon);
   
   procedure Snap_To_Crystallographic_Translations
     (M : in out Symmetry_Operator) is
   begin
      for I in 1 .. M'Last(1) - 1 loop
         M (I,4) := M (I,4) - Float'Floor (M (I,4));
         if abs (M (I,4) - 1.0/3.0) < Eps  then
            M (I,4) := 1.0/3.0;
         elsif abs (M (I,4) - 2.0/3.0) < Eps then
            M (I,4) := 2.0/3.0;
         elsif abs (M (I,4) - 1.0/6.0) < Eps then
            M (I,4) := 1.0/6.0;
         elsif abs (M (I,4) - 5.0/6.0) < Eps then
            M (I,4) := 5.0/6.0;
         elsif abs (M (I,4) - 1.0/2.0) < Eps then
            M (I,4) := 1.0/2.0;
         elsif abs (M (I,4) - 1.0/4.0) < Eps then
            M (I,4) := 1.0/4.0;
         elsif abs (M (I,4) - 3.0/4.0) < Eps then
            M (I,4) := 3.0/4.0;
         elsif abs (M (I,4) - 0.0) < Eps or else abs (M (I,4) - 1.0) < Eps then
            M (I,4) := 0.0;
         end if;
      end loop;
   end Snap_To_Crystallographic_Translations;
   
   procedure Add
     (M : in out Symmetry_Operator; T : Crystallographic_Translation) is
   begin
      for I in 1..3 loop
         M (I,4) := M (I,4) + 
           Float (T (I).Numerator) / Float (T (I).Denominator);
      end loop;
      Snap_To_Crystallographic_Translations (M);
   end;
   
   type Matrix3x3 is array (1..3,1..3) of Float;
   
   function Det (M : Matrix3x3) return Float is
   ( 
     M(1,1)*M(2,2)*M(3,3) + 
       M(1,2)*M(2,3)*M(3,1) + 
       M(2,1)*M(3,2)*M(1,3) -
       M(1,3)*M(2,2)*M(3,1) - 
       M(2,1)*M(1,2)*M(3,3) - 
       M(1,1)*M(2,3)*M(3,2)
   );
      
   function Invert (M : Matrix3x3) return Matrix3x3 is
      D : Float := Det (M);
      
      function Adjunct (P, Q : Integer) return Float is
         A : array (1..4,1..2) of Float;
         K, L : Integer;
         Coeff : Float;
      begin
         K := 1; L := 1;
         for I in M'Range(1) loop
            if I /= P then 
               for J in M'Range(2) loop
                  if J /= Q then
                     A(K,L) := M(I,J);
                     L := L + 1;
                  end if;
               end loop;
               L := 1;
               K := K + 1;
            end if;
         end loop;
         if (P + Q) mod 2 = 0 then
            Coeff := 1.0;
         else
            Coeff := -1.0;
         end if;
         return Coeff * (A(1,1)*A(2,2) - A(1,2)*A(2,1));
      end Adjunct;
      
   begin -- Invert
      return (
              ( Adjunct(1,1)/D, Adjunct(2,1)/D, Adjunct(3,1)/D ),
              ( Adjunct(1,2)/D, Adjunct(2,2)/D, Adjunct(3,2)/D ),
              ( Adjunct(1,3)/D, Adjunct(2,3)/D, Adjunct(3,3)/D )
             );
   end Invert;
   
   -- a very specific inversion routine for sympos:
   function Invert (S : Symmetry_Operator) return Symmetry_Operator is
      R : Matrix3x3;
      Inv : Symmetry_Operator;
   begin
      for I in R'Range(1) loop
         for J in R'Range(2) loop
            R (I,J) := S (I,J);
         end loop;
      end loop;
      
      -- rotation matrix of the inverse symop is an inverted rotation
      --  matrix:
      R := Invert (R);
      
      for I in R'Range(1) loop
         for J in R'Range(2) loop
            Inv (I,J) := R (I,J);
         end loop;
      end loop;
      
      -- compute the inverse translation:
      
      -- assume R' is the inverse of R:
      -- R * R' = I, where I is the unity matrix.
      -- then:
      -- (R,t) * (R',t') = (R*R', R't + t') = (I,0)
      -- =>
      -- R't + t' = 0
      -- =>
      -- t' = -R' * t
      
      for I in R'Range(1) loop
         Inv (I,4) := 0.0;
         for J in R'Range(2) loop
            Inv (I,4) := Inv (I,4) - R (I,J) * S (J,4);
         end loop;
      end loop;
      
      -- The last row is (0,0,0,1):
      Inv (4,1) := 0.0;
      Inv (4,2) := 0.0;
      Inv (4,3) := 0.0;
      Inv (4,4) := 1.0;
      
      return Inv;
   end Invert;
   
   type Float_Matrix is array (Integer range <>, Integer range <>) of Float;
   
   function Transpose (A : Symmetry_Operator) return Symmetry_Operator is
      R : Symmetry_Operator;
   begin
      -- Transpose the rotation part:
      for I in 1 .. 3 loop
         for J in 1 .. 3 loop
            R (I, J) := A(J, I);
         end loop;
      end loop;
      -- Copy the rest:
      for I in 1 .. 3 loop
         R(I, 4) := A (I, 4);
         R(4, I) := A (4, I);
      end loop;
      -- The (4, 4) elemet is always the same:
      R (4, 4) := 1.0;
      return R;
   end Transpose;
   
   function "*" (M1, M2 : Symmetry_Operator) return Symmetry_Operator is
      M : Symmetry_Operator;
   begin
      pragma Assert (M1'Last(2) = M2'Last(1));
      
      for I in M1'Range(1) loop
         for J in M2'Range(2) loop
            M (I,J) := 0.0;
            for K in M2'Range(1) loop
               M (I,J) := M (I,J) + M1 (I,K) * M2(K,J);
            end loop;
         end loop;
      end loop;
      Snap_To_Crystallographic_Translations (M);
      return M;
   end "*";
   
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
   
   function Has_Symmetry_Operator
     (
      Symmetry_Operators : Symmetry_Operator_Array;
      Last_Symmetry_Operator_Index : Positive;
      Lookup_Symmetry_Operator : Symmetry_Operator
     )
     return Boolean is
   begin
      for I in 1 .. Last_Symmetry_Operator_Index loop
         if Symmetry_Operators (I) = Lookup_Symmetry_Operator then
            return True;
         end if;
      end loop;
      return False;
   end;
   
   procedure Apply_Change_Of_Basis
     (
      Symmetry_Operators : in out Symmetry_Operator_Array;
      N_Symmetry_Operators : in out Positive;
      Change_Of_Basis : in Symmetry_Operator;
      Debug_Print_Matrices : Boolean := False
     ) is
      
      Max_Centering : constant Positive := 8;
      Centering : Symmetry_Operator_Array (1..Max_Centering);
      N_Centering : Natural := 0;
      
   begin
      if Change_Of_Basis /= Unity_Matrix then
         declare
            V : Symmetry_Operator := Change_Of_Basis;
            V_Inv : Symmetry_Operator;
            New_Symmetry_Operators : Symmetry_Operator_Array (1 .. N_Symmetry_Operators);
            N_New_Operators : Positive := 1;
            New_Operator : Symmetry_Operator;
         begin
            New_Symmetry_Operators (1) := Symmetry_Operators (1);
            V_Inv := Invert (V);
            for I in 2..N_Symmetry_Operators loop
               New_Operator := V * Symmetry_Operators (I) * V_Inv;
               if not Has_Symmetry_Operator
                 (
                  New_Symmetry_Operators,
                  N_New_Operators,
                  New_Operator
                 ) then
                  N_New_Operators := N_New_Operators + 1;
                  New_Symmetry_Operators (N_New_Operators) := New_Operator;
               end if;
            end loop;
            N_Symmetry_Operators := N_New_Operators;
            Symmetry_Operators (2 .. N_New_Operators) :=
              New_Symmetry_Operators (2 .. N_New_Operators);
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
                  R (I) := R (I) + S (I,K) * T (K);
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
            
      begin -- declare
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
         
      -- Add centering matrices:
         
      declare
         M : Natural := N_Symmetry_Operators;
         New_Symmetry_Operator : Symmetry_Operator;
      begin
         for C in 1..N_Centering loop
            for S in 1..N_Symmetry_Operators loop
               New_Symmetry_Operator :=
                 Symmetry_Operators (S) * Centering (C);
               if not Has_Symmetry_Operator (Symmetry_Operators, M, 
                                             New_Symmetry_Operator) then
                  M := M + 1;
                  Symmetry_Operators (M) := New_Symmetry_Operator;
               end if;
            end loop;
         end loop;
         N_Symmetry_Operators := M;
      end;
      
      -- Print out matrices for debug purposes if requested:
      
      -- Print out all matrices if requested:
      
      if Debug_Print_Matrices then
         Put_Line (Standard_Error, "Rotations:");
         for I in 1..N_Symmetry_Operators loop
            Put_Line (Standard_Error, Integer'Image (I));
            Put (Standard_Error, Symmetry_Operators (I));
            New_Line (Standard_Error);
         end loop;
         
         Put_Line (Standard_Error, "Centerings:");
         for I in 1..N_Centering loop
            Put (Standard_Error, Centering (I));
            New_Line (Standard_Error);
         end loop;
         New_Line (Standard_Error);
         
         Put_Line (Standard_Error, "Change of basis:");
         Put (Standard_Error, Change_Of_Basis);
         
         Put_Line (Standard_Error, "Inverse Change of basis:");
         Put (Standard_Error, Invert (Change_Of_Basis));
         New_Line (Standard_Error);
      end if;
      
   end Apply_Change_Of_Basis;
   
   procedure Build_Group
     (
      Operators : in out Symmetry_Operator_Array;
      N_Operators : in out Natural
     )
   is
      N, M : Natural := N_Operators;
      New_Operator : Symmetry_Operator;
   begin
      loop
         for I in 2 .. N loop
            New_Operator := Operators (I) * Operators (N);
            if not Has_Symmetry_Operator (Operators, M, New_Operator) then
               M := M + 1;
               if M > Operators'Last then
                  for I in Operators'Range loop
                     Put (As_String (Operators (I)));
                     New_Line;
                     Flush;
                  end loop;
                  raise CONSTRAINT_ERROR
                    with "Operator no." & M'Image & " exceeds capacity " &
                    "of the array when storing " &  As_String (New_Operator);
               end if;
               Operators (M) := New_Operator;
            end if;
         end loop;
         N := N + 1;
         exit when N > M or else M >= Operators'Last;
      end loop;
      N_Operators := M;
   end Build_Group;

end Symmetry_Operations;
