with Ada.Text_IO; use Ada.Text_IO;

with Change_Of_Basis;           use Change_Of_Basis;
with HM_Symbols;                use HM_Symbols;
with Shmueli_Symbol_Parser;     use Shmueli_Symbol_Parser;

package body HM_Symbol_Parser is
   
   -- The (a,b,c) change-of-basis notation is taken from [1], p 100.
   
   -- [1] Zwart, P. H.; Grosse-Kunstleve, R. W.; Lebedev, A. A.;
   --  Murshudov, G. N. & Adams, P. D. (2007) Surprises and pitfalls
   --  arising from (pseudo)symmetry. Acta Crystallographica Section D
   --  Biological Crystallography 64(1), 99-107. International Union
   --  of Crystallography (IUCr). DOI:
   --  https://doi.org/10.1107/s090744490705531x
   
   procedure Inc (X : in out Positive) is
   begin
      X := X + 1;
   end;
   
   pragma Inline (Inc);
   
   procedure Lookup_H_M_Symbol_Rotations
     (
      Symbol : in String;
      Pos : in out Positive;
      Rotations : out Symmetry_Operator_Array;
      N_Rotations : in out Natural
     )
   is
      Start : Positive := Pos;
      J : Positive := 1;
   begin
      N_Rotations := 0;

      while Pos <= Symbol'Length and then Symbol (Pos) /= '(' loop
         Inc (Pos);
      end loop;
      
      declare
         HM_Part : String (1 .. Pos - Start);
         HM_Symbol : HM_Symbol_Type;
      begin
         for I in Start .. Pos - 1 loop
            if Symbol (I) /= ' ' then
               HM_Part (J) := Symbol (I);
               Inc (J);
            end if;
         end loop;
         
         HM_Symbol := +HM_Part (1 .. J - 1);
         
         for I in HM_Symbol_Table'Range loop
            if HM_Symbol = HM_Symbol_Table(I).HM_Symbol then
               declare
                  Matrices : Symmetry_Operator_Array :=
                    Decode_Shmueli_Symbol (String (HM_Symbol_Table(I).Shmueli_Symbol));
               begin      
                  N_Rotations := Matrices'Length;
                  Rotations (1 .. N_Rotations) := Matrices;
               end;
               return;
            end if;
         end loop;
         
         raise HM_SYMBOL_NOT_FOUND with
           "Hermann-Mauguin symbol '" & Symbol & "' " &
           "was not found in internal tables";
      end;
   end;
   
   ----------------------------------------------------------------------------
   
   function As_String (S : Symmetry_Operator) return String;

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
   
   function Decode_Hermann_Mauguin_Symbol (Symbol : in String) 
                                          return Symmetry_Operator_Array is
      
      Max_Symmetry_Operators : constant Integer := 192;
      
      Symmetry_Operators :
        Symmetry_Operator_Array (1 .. Max_Symmetry_Operators);
      N_Symmetry_Operators : Natural := 0;
      
      Pos : Positive := 1; -- current position in the string 'Symbol'.
      
      Change_Of_Basis : Symmetry_Operator;
      
      Max_Centering : constant Positive := 8;
      Centering : Symmetry_Operator_Array (1..Max_Centering);
      N_Centering : Natural := 0;
      
   begin
      
      Lookup_H_M_Symbol_Rotations (Symbol, Pos, Symmetry_Operators,
                                   N_Symmetry_Operators);
      
      Get_Change_Of_Basis (Symbol, Pos, Change_Of_Basis);
      
      -- Apply the change-of-basis operator:
      
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
      
      -- Reconstruct all rotation operators once again:
      
      Build_Group (Symmetry_Operators, N_Symmetry_Operators);
      
      return Symmetry_Operators (1..N_Symmetry_Operators);
   end Decode_Hermann_Mauguin_Symbol;
   
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
   
end HM_Symbol_Parser;
