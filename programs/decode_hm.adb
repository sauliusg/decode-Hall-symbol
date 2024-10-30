pragma Ada_2022;

with Text_IO;                   use Text_IO;
with Ada.Integer_Text_IO;       use Ada.Integer_Text_IO;
with Ada.Command_Line;          use Ada.Command_Line;
with Ada.Environment_Variables; use Ada.Environment_Variables;
with Ada.Strings.Fixed;         use Ada.Strings.Fixed;
with Ada.Strings.Maps;          use Ada.Strings.Maps;

with Symmetry_Operations;       use Symmetry_Operations;
with HM_Symbols;                use HM_Symbols;
with Shmueli_Symbol_Parser;     use Shmueli_Symbol_Parser;

with Project_Version;

procedure Decode_HM is
   
   -- This program decodes Hermann-Mauguin Crystallographic space
   --  group symbols and outputs space group symmetry operators as
   --  lists of general position coordinates.
   --
   
   -- The (a,b,c) change-of-basis notation is taken from [1], p 100.
   
   -- [1] Zwart, P. H.; Grosse-Kunstleve, R. W.; Lebedev, A. A.;
   --  Murshudov, G. N. & Adams, P. D. (2007) Surprises and pitfalls
   --  arising from (pseudo)symmetry. Acta Crystallographica Section D
   --  Biological Crystallography 64(1), 99-107. International Union
   --  of Crystallography (IUCr). DOI:
   --  https://doi.org/10.1107/s090744490705531x
   
   Debug_Print_Matrices : Boolean := False;
   
   function IDENTITY return Axis_Order_Type 
     renames Symmetry_Operations.IDENTITY;
   
   function Machine_Epsilon return Float is
      Epsilon : Float := 1.0;
   begin
      while 1.0 + Epsilon /= 1.0 loop
         Epsilon := Epsilon / 2.0;
      end loop;
      return Epsilon;
   end;
   
   Eps : constant Float := 16.0 * Machine_Epsilon;

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
   
   -- -------------------------------------------------------------------------
   
   procedure Inc (X : in out Positive) is
   begin
      X := X + 1;
   end;
   
   pragma Inline (Inc);
   
   HM_SYMBOL_NOT_FOUND : exception;
   
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
   
   
   subtype Character_Set is Ada.Strings.Maps.Character_Set;
   
   procedure Expect (
                     Symbol : in String;
                     Pos : in out Integer;
                     Ch_Set : in Character_Set
                    ) 
   is
   begin
      Skip_Spaces (Symbol, Pos);
      if Pos <= Symbol'Last then
         if not Is_In( Symbol (Pos), Ch_Set) then
            raise UNEXPECTED_SYMBOL with
              "symbol " & Character'Image (Symbol (Pos)) & " " &
              "is not expected at position" & Integer'Image (Pos) &
              " in """ & Symbol & """" &
              ", expecting one of """ &
              To_Sequence (Ch_Set) & """";
         end if;
      else
         raise UNEXPECTED_SYMBOL with
           "unexpected end-of-string";
      end if;
   end Expect;
   
   procedure Skip (
                   Symbol : in String;
                   Pos : in out Integer;
                   Ch_Set : in Character_Set
                  )
   is
   begin
      Expect (Symbol, Pos, Ch_Set);
      while Pos <= Symbol'Length and then Is_In (Symbol (Pos), Ch_Set) loop
         Pos := Pos + 1;
      end loop;
   end;

   ----------------------------------------------------------------------------
   -- A simple recursive descent parser for the change-of-basis operator:
   
   -- parse a fractional number (e.g. "2/3") and return its value as Float:
   procedure Inc (Result : in out Float; D : in Float) is
   begin
      Result := Result + D;
   end;
   
   function Get_Number (Symbol : in String; Pos : in out Integer) return Float
   is
      Numerator : Natural := 1;
      Denominator : Natural := 1;
      Final_Pos : Integer := Pos;
   begin
      while Final_Pos <= Symbol'Last and then 
        Symbol (Final_Pos) in '0'..'9'
      loop
         Final_Pos := Final_Pos + 1;
      end loop;
      Numerator := Integer'Value (Symbol (Pos..Final_Pos-1));
      Skip_Spaces (Symbol, Final_Pos);
      Pos := Final_Pos;
      if Pos <= Symbol'Last and then Symbol (Pos) = '/' then
         Pos := Pos + 1;
         Skip_Spaces (Symbol, Pos);
         Final_Pos := Pos;
         while Final_Pos <= Symbol'Last and then 
           Symbol (Final_Pos) in '0'..'9'
         loop
            Final_Pos := Final_Pos + 1;
         end loop;
         Denominator := Integer'Value (Symbol (Pos..Final_Pos-1));
         Pos := Final_Pos;
      end if;
      return Float (Numerator) / Float (Denominator);
   end Get_Number;
   
   -- parse the '+2*a' factor:
   procedure Parse_Factor
     (
      Symbol : in String;
      Pos : in out Integer;
      Change_Of_Basis : out Symmetry_Operator;
      Row : in Integer;
      Factor : in Float
     ) is
   begin
      Skip_Spaces (Symbol, Pos);
      Expect (Symbol, Pos, To_Set ("abc"));
      case Symbol (Pos) is
         when 'a' => Change_Of_Basis (Row, 1) := Factor;
         when 'b' => Change_Of_Basis (Row, 2) := Factor;
         when 'c' => Change_Of_Basis (Row, 3) := Factor;
         when others =>
            raise UNEXPECTED_SYMBOL with
              "unexpected character " & Character'Image (Symbol (Pos));
      end case;
      Pos := Pos + 1;
   end Parse_Factor;
   
   -- Parse the "+a", "b", "1/4" parts in the "+a+2*b+1/4:
   procedure Parse_Term
     (
      Symbol : in String;
      Pos : in out Integer;
      Change_Of_Basis : out Symmetry_Operator;
      Row : in Integer
     ) is
      Factor : Float := 1.0;
   begin
      Skip_Spaces (Symbol, Pos);
      
      if Pos <= Symbol'Last then
         if Symbol (Pos) = '+' then
            Pos := Pos + 1;
         elsif Symbol (Pos) = '-' then
            Factor := -1.0;
            Pos := Pos + 1;
         end if;
      end if;
      
      Expect (Symbol, Pos, To_Set ("0123456789abc"));
      
      if Pos <= Symbol'Last then
         case Symbol (Pos) is
            when 'a' => 
               Change_Of_Basis (Row, 1) := Factor;
               Pos := Pos + 1;
            when 'b' => 
               Change_Of_Basis (Row, 2) := Factor;
               Pos := Pos + 1;
            when 'c' => 
               Change_Of_Basis (Row, 3) := Factor;
               Pos := Pos + 1;
            when '0'..'9' =>
               Factor := Factor * Get_Number (Symbol, Pos);
               Skip_Spaces (Symbol, Pos);
               if Pos <= Symbol'Length and then Symbol (Pos) = '*' then
                  Pos := Pos + 1;
                  Parse_Factor (Symbol, Pos, Change_Of_Basis, Row, Factor);
               else
                  Inc (Change_Of_Basis (Row, 4), Factor);
               end if;
            when others =>
               raise UNEXPECTED_SYMBOL with
                 "unexpected symbol " & Character'Image (Symbol (Pos)) &
                 " in the symop """ & Symbol & """";
         end case;
      end if;
   end Parse_Term;
   
   -- parse the "2*b+1/4" part in the "2*b+1/4,c,a-1/3" operator:
   procedure Parse_Symmetry_Operator_Component
     (
      Symbol : in String;
      Pos : in out Integer;
      Change_Of_Basis : out Symmetry_Operator;
      Row : Integer
     ) is
   begin
      loop
         Parse_Term (Symbol, Pos, Change_Of_Basis, Row);
         Skip_Spaces (Symbol, Pos);
         if Pos > Symbol'Length or else 
           Is_In (Symbol (Pos), To_Set (",)")) then
            exit;
         end if;
      end loop;
   end Parse_Symmetry_Operator_Component;
   
   procedure Interpret_Change_Of_Basis_Matrix
      (
       Symbol : in String;
       Pos : in out Integer;
       Change_Of_Basis : out Symmetry_Operator
      )
   is
   begin
      Change_Of_Basis := Zero_Matrix;
      Change_Of_Basis (4,4) := 1.0;
      Skip (Symbol, Pos, To_Set('('));
      Parse_Symmetry_Operator_Component (Symbol, Pos, Change_Of_Basis, 1);
      Skip (Symbol, Pos, To_Set(','));
      Parse_Symmetry_Operator_Component (Symbol, Pos, Change_Of_Basis, 2);
      Skip (Symbol, Pos, To_Set(','));
      Parse_Symmetry_Operator_Component (Symbol, Pos, Change_Of_Basis, 3);
      Skip (Symbol, Pos, To_Set(')'));
   end Interpret_Change_Of_Basis_Matrix;
   
   procedure Get_Shift_Of_Origin (
                                  Symbol : in String;
                                  Pos : in out Integer;
                                  Change_Of_Basis : out Symmetry_Operator
                                 )
   is
      Shift : Integer;
      Sign : Integer := 1;
      S : Symmetry_Operator := Unity_Matrix;
   begin
      Skip_Spaces (Symbol, Pos);
      
      if Pos <= Symbol'Last and then Symbol (Pos) = '(' then
         Pos := Pos + 1;
         
         for I in 1 .. 3 loop
            Expect (Symbol, Pos, To_Set ("-0123456789"));
            if Symbol (Pos) = '-' then
               Sign := -1;
               Pos := Pos + 1;
            end if;
            Get (Symbol (Pos..Symbol'Last), Shift, Pos);
            Pos := Pos + 1;
            S (I,4) := Float (Sign * Shift) / 12.0;
         end loop;
         
         Expect (Symbol, Pos, To_Set (')'));         
      end if;
      Change_Of_Basis := S;
   end Get_Shift_Of_Origin;
   
   function Has_Only_Characters (S : String; CS : Character_Set) return Boolean
   is
   begin
      for C of S loop
         if not Is_In (C, CS) then
            return False;
         end if;
      end loop;
      return True;
   end;
   
   procedure Get_Change_Of_Basis (
                                  Symbol : in String;
                                  Pos : in out Integer;
                                  Change_Of_Basis : out Symmetry_Operator
                                 )
   is
   begin
      if Has_Only_Characters (Symbol (Pos..Symbol'Last),
                              To_Set ("-( 0123456789)")) then
         Get_Shift_Of_Origin (Symbol, Pos, Change_Of_Basis);
      else
         Interpret_Change_Of_Basis_Matrix (Symbol, Pos, Change_Of_Basis);
      end if;
   end;
   
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
   end;
   
   function Decode_Hermann_Mauguin (Symbol : in String) 
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
            V : Symmetry_Operator := Invert (Transpose (Change_Of_Basis));
            V_Inv : Symmetry_Operator;
         begin
            V_Inv := Invert (V);
            for I in 2..N_Symmetry_Operators loop
               Symmetry_Operators (I) := V * Symmetry_Operators (I) * V_Inv;
            end loop;
         end;
      end if;
      
      -- Generate additional centering operators:
      
      declare
         type Vector_Components is array (1..4) of Float;
         type Vector_Type is record
            Value : Vector_Components;
         end record;
         
         function "*" (S : Symmetry_Operator; T : Vector_Type)
                      return Vector_Type 
         is
            R : Vector_Type := (Value => (others => 0.0));
         begin
            for I in R.Value'Range loop
               for K in T.Value'Range loop
                  R.Value (I) := R.Value (I) +
                    S (I,K) * T.Value (K);
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
               S (I,4) := T.Value (I);
            end loop;
            Snap_To_Crystallographic_Translations (S);
            return S;
         end;
         
         function Is_Centering (T : Vector_Type) return Boolean is            
            function Fract (X : Float) return Float is (X - Float'Floor (X));
         begin
            for Component of T.Value loop
               if abs (Fract (Component)) >= Eps then
                  return True;
               end if;
            end loop;
            return False;
         end;
         
         Unit_Vectors : array (1..3) of Vector_Type :=
           (
            (Value => (1.0, 0.0, 0.0, 1.0)),
            (Value => (0.0, 1.0, 0.0, 1.0)),
            (Value => (0.0, 0.0, 1.0, 1.0))
           );
         
         C_O_B_Rotation : Symmetry_Operator := Invert (Change_Of_Basis);
         
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
         Put (Standard_Error, Transpose (Change_Of_Basis));
         Put_Line (Standard_Error, "Inverse Change of basis:");
         Put (Standard_Error, Invert (Transpose (Change_Of_Basis)));
         New_Line (Standard_Error);
      end if;
      
      -- Reconstruct all rotation operators once again:
      
      Build_Group (Symmetry_Operators, N_Symmetry_Operators);
      
      return Symmetry_Operators (1..N_Symmetry_Operators);
   end Decode_Hermann_Mauguin;
   
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
   
begin
   
   if Exists ("DECODE_HM_DEBUG") and then
     (
      Value ("DECODE_HM_DEBUG") = "1" or else
        Value ("DECODE_HM_DEBUG") = "true"
     )
   then
      Debug_Print_Matrices := True;
   end if;
   
   for I in 1 .. Argument_Count loop
      if I > 1 then
         New_Line;
      end if;
            
      if Index ("--help", Argument (I)) = 1
      then
         Put_Line ("Decode a Hermann-Mauguin Crystallographic spacegroup symbol " &
                     "given on the command line");
         Put_Line ("and output symmtetry operators as general position " &
                     "point coordinates, e.g. '-X,-Y,Z+1/2'");
         New_Line;
         Put_Line ("USAGE:");
         Put_Line ("  " & Command_Name & " 'P 21 21 21'");
         Put_Line ("  " & Command_Name & " --help");
      elsif Index ("--version", Argument (I)) = 1 then
         Put (Command_Name & " " & Project_Version.Version);
         if Project_Version.VCS_Text /= "" then
            Put (" " & Project_Version.VCS_Text);
         end if;
         New_Line;
      else
         if Debug_Print_Matrices then
            Put_Line (Standard_Error, Argument (I));
         end if;
         
         declare
            Symmetry_Operators : Symmetry_Operator_Array :=
              Decode_Hermann_Mauguin (Argument (I));
         begin
            if Debug_Print_Matrices then
               Put_Line (Standard_Error, "Symmetry_Operators:");
               for I in Symmetry_Operators'Range loop
                  Put (Standard_Error, Symmetry_Operators (I));
                  New_Line (Standard_Error);
               end loop;
            end if;
            
            for I in Symmetry_Operators'Range loop
               Put_Line (As_String (Symmetry_Operators (I)));
            end loop;
         end;
      end if;
   end loop;

end Decode_HM;
