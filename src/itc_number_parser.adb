with Ada.Text_IO;           use Ada.Text_IO;
with Ada.Integer_Text_IO;   use Ada.Integer_Text_IO;
with Ada.Strings.Fixed;     use Ada.Strings.Fixed;

with ITC_Number_Tables;     use ITC_Number_Tables;
with Shmueli_Symbol_Parser; use Shmueli_Symbol_Parser;

package body ITC_Number_Parser is   

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
      
   begin
      Symmetry_Operators (1) := Unity_Matrix;
      
      if Index (":", ITC_Number_And_Cell) = 0 then
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
      end if;
      
      -- Reconstruct all symmetry operators:
      
      Build_Group (Symmetry_Operators, N_Symmetry_Operators);
      
      return Symmetry_Operators (1 .. N_Symmetry_Operators);
   end;
   
end ITC_Number_Parser;
