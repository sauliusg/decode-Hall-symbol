with Ada.Text_IO; use Ada.Text_IO; 
package body Parser_Tools is 
   
   procedure Skip_Spaces (S : in String; Pos : in out Integer ) is
   begin
      while Pos <= S'Last and then S (Pos) = ' ' loop
         Pos := Pos + 1;
      end loop;
   end;
   
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
   
   function Has_Characters (S : String; CS : Character_Set) return Boolean
   is
   begin
      for C of S loop
         if Is_In (C, CS) then
            return True;
         end if;
      end loop;
      return False;
   end;

end Parser_Tools;
