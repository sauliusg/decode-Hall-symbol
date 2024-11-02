with Ada.Text_IO; use Ada.Text_IO;

with Change_Of_Basis;           use Change_Of_Basis;
with HM_Tables;                 use HM_Tables;
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
   
   function Decode_Hermann_Mauguin_Symbol (Symbol : in String) 
                                          return Symmetry_Operator_Array is
      
      Max_Symmetry_Operators : constant Integer := 192;
      
      Symmetry_Operators :
        Symmetry_Operator_Array (1 .. Max_Symmetry_Operators);
      N_Symmetry_Operators : Natural := 0;
      
      Pos : Positive := 1; -- current position in the string 'Symbol'.
      
      Change_Of_Basis : Symmetry_Operator;
      
   begin
      
      Lookup_H_M_Symbol_Rotations (Symbol, Pos, Symmetry_Operators,
                                   N_Symmetry_Operators);
      
      Get_Change_Of_Basis (Symbol, Pos, Change_Of_Basis);
      
      -- Apply the change-of-basis operator:
      
      Apply_Change_Of_Basis
        (
         Symmetry_Operators,
         N_Symmetry_Operators,
         Change_Of_Basis,
         Debug_Print_Matrices
        );
      
      -- Reconstruct all rotation operators once again:
      
      Build_Group (Symmetry_Operators, N_Symmetry_Operators);
      
      return Symmetry_Operators (1..N_Symmetry_Operators);
   end Decode_Hermann_Mauguin_Symbol;
   
end HM_Symbol_Parser;
