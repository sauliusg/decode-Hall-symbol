with Ada.Integer_Text_IO;  use Ada.Integer_Text_IO;
with Ada.Strings.Fixed;    use Ada.Strings.Fixed;
with Ada.Strings.Maps;     use Ada.Strings.Maps;

with Parser_Tools;         use Parser_Tools;

package body Change_Of_Basis is
   
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
      Expect (Symbol, Pos, To_Set ("abcxyz"));
      case Symbol (Pos) is
         when 'a' | 'x' => Change_Of_Basis (Row, 1) := Factor;
         when 'b' | 'y' => Change_Of_Basis (Row, 2) := Factor;
         when 'c' | 'z' => Change_Of_Basis (Row, 3) := Factor;
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
      
      Expect (Symbol, Pos, To_Set ("0123456789abcxyz"));
      
      if Pos <= Symbol'Last then
         case Symbol (Pos) is
            when 'a' | 'x' => 
               Change_Of_Basis (Row, 1) := Factor;
               Pos := Pos + 1;
            when 'b' | 'y' => 
               Change_Of_Basis (Row, 2) := Factor;
               Pos := Pos + 1;
            when 'c' | 'z' => 
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
   
   procedure Get_Change_Of_Basis (
                                  Symbol : in String;
                                  Pos : in out Integer;
                                  Change_Of_Basis : out Symmetry_Operator
                                 )
   is
      Is_Inverted_Matrix, Is_Direct_Matrix : Boolean;
   begin
      if Has_Only_Characters (Symbol (Pos..Symbol'Last),
                              To_Set ("-( 0123456789)")) then
         Get_Shift_Of_Origin (Symbol, Pos, Change_Of_Basis);
      else
         Is_Inverted_Matrix :=
           Has_Characters (Symbol (Pos..Symbol'Last), To_Set ("abc"));
         
         Is_Direct_Matrix :=
           Has_Characters (Symbol (Pos..Symbol'Last), To_Set ("xyz"));
         
         if Is_Inverted_Matrix and Is_Direct_Matrix then
            raise WRONG_CHANGE_OF_BASIS with
              "the change of basis operator may contain either only " &
              "the ""abc"" characters or only the ""xyz"" characters, " &
              "but """ & Symbol (Pos..Symbol'Last) & """ contains both";
         end if;
         
         Interpret_Change_Of_Basis_Matrix (Symbol, Pos, Change_Of_Basis);
         
         if Is_Inverted_Matrix then
            Change_Of_Basis := Invert (Transpose (Change_Of_Basis));
         end if;
      end if;
   end;
   
end Change_Of_Basis;
