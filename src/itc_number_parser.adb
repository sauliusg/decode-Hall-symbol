with Ada.Text_IO;           use Ada.Text_IO;
with Ada.Integer_Text_IO;   use Ada.Integer_Text_IO;
with Ada.Strings.Fixed;     use Ada.Strings.Fixed;

with Change_Of_Basis;       use Change_Of_Basis;
with ITC_Number_Tables;     use ITC_Number_Tables;
with Shmueli_Symbol_Parser; use Shmueli_Symbol_Parser;

package body ITC_Number_Parser is   
   
   function Lookup_ITC_Number_And_Setting (S : String; Pos : in out Positive)
                                          return Shmueli_Symbol_Type is
      ITC_NC : ITC_Number_And_Cell_Symbol_Type;
      Start : Positive := Pos;      
   begin
      while Pos <= S'Last and then S (Pos) /= ' ' loop
         Pos := Pos + 1;
      end loop;
      
      ITC_NC := +S(Start .. Pos - 1);
      
      for I in ITC_Number_And_Cell_Symbol_Table'Range loop
         if ITC_NC = 
           ITC_Number_And_Cell_Symbol_Table (I).
           ITC_Number_And_Cell_Symbol then
            return ITC_Number_And_Cell_Symbol_Table (I).Shmueli_Symbol;
         end if;
      end loop;
      
      raise ITC_NUMBER_SYMBOL_NOT_FOUND with
        "ITC symbol """ & S(Start .. Pos - 1) & """ could not be " &
        "found in the lookup tables";
   end;
   
   function Lookup_ITC_Number (SG_Number : in Space_Group_Number_Type)
                              return Symmetry_Operator_Array is
   begin
      return Decode_Shmueli_Symbol
        (String (ITC_Number_Shmueli_Symbols (SG_Number)));
   end;
   
   function Parse_ITC_Number (ITC_Number_And_Cell : in String)
                              return Symmetry_Operator_Array
   is
      -- current position in the string 'ITC_Number_And_Cell':
      Pos : Positive := 1;
      
      Max_Symmetry_Operators : constant Integer := 192;
      
      Symmetry_Operators :
        Symmetry_Operator_Array (1 .. Max_Symmetry_Operators);
      
      N_Symmetry_Operators : Positive := 1;
      
      Change_Of_Basis : Symmetry_Operator;
      
      Semicolon_Index : Natural := Index (ITC_Number_And_Cell, ":");
   begin
      Symmetry_Operators (1) := Unity_Matrix;
      
      if Semicolon_Index = 0 then
         declare
            SG_Number : Space_Group_Number_Type;
         begin
            Get (ITC_Number_And_Cell, SG_Number, Pos);
            declare
               Default_Symmetry_Operators : Symmetry_Operator_Array :=
                 Lookup_ITC_Number (SG_Number);
               N_Default_Operators : Positive := Default_Symmetry_Operators'Last;
            begin
               Symmetry_Operators (1 .. N_Default_Operators) :=
                 Default_Symmetry_Operators;
               N_Symmetry_Operators := N_Default_Operators;
            end;
         end;
      else
         if Semicolon_Index < ITC_Number_And_Cell'Length and then
           ITC_Number_And_Cell (Semicolon_Index + 1) = '1' and then
           (Semicolon_Index = ITC_Number_And_Cell'Length - 1 or else
              not (ITC_Number_And_Cell (Semicolon_Index + 2) in '0' .. '9'))
         then
            declare
               SG_Number : Integer;
               Start : Integer := ITC_Number_And_Cell'First;
            begin
               Get (ITC_Number_And_Cell (Start .. Semicolon_Index - 1), 
                    SG_Number, Pos);
               declare
                  Default_Symmetry_Operators : Symmetry_Operator_Array :=
                    Lookup_ITC_Number (SG_Number);
                  N_Default_Operators : Positive := 
                    Default_Symmetry_Operators'Last;
               begin
                  Symmetry_Operators (1 .. N_Default_Operators) :=
                    Default_Symmetry_Operators;
                  N_Symmetry_Operators := N_Default_Operators;
               end;
            end;
         else
            declare
               SH : Shmueli_Symbol_Type :=
                 Lookup_ITC_Number_And_Setting (ITC_Number_And_Cell, Pos);
               
               Matrices : Symmetry_Operator_Array :=
                 Decode_Shmueli_Symbol (String (SH));
            begin
               N_Symmetry_Operators := Matrices'Length;
               Symmetry_Operators (1 .. Matrices'Length) := Matrices;
            end;
         end if;
         
         Pos := Semicolon_Index;
         while Pos <= ITC_Number_And_Cell'Last and then
           ITC_Number_And_Cell (Pos) /= ' ' loop
            Pos := Pos + 1;
         end loop;
         
         Get_Change_Of_Basis (ITC_Number_And_Cell, Pos, Change_Of_Basis);
         
         Apply_Change_Of_Basis
           (
            Symmetry_Operators,
            N_Symmetry_Operators,
            Change_Of_Basis,
            Debug_Print_Matrices
           );
      end if;
      
      
      -- Reconstruct all symmetry operators:
      
      Build_Group (Symmetry_Operators, N_Symmetry_Operators);
      
      return Symmetry_Operators (1 .. N_Symmetry_Operators);
   end;
   
end ITC_Number_Parser;
