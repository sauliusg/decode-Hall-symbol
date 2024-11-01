with Ada.Strings.Maps;          use Ada.Strings.Maps;

package Parser_Tools is 
   
   -- A set of simple routines to help with constructing a simple
   --  string parser.
   
   subtype Character_Set is Ada.Strings.Maps.Character_Set;
   
   UNEXPECTED_SYMBOL : exception;
   
   -- Advance the 'Pos' index in 'S' past the space characters:
   
   procedure Skip_Spaces (S : in String; Pos : in out Integer );

   -- Check that the next character of the 'Symbol' is in the
   --  specified set 'Ch_Set' and advance the position index
   --  'Pos'. Raise an UNEXPECTED_SYMBOL exception if the character
   --  was not an expected one:
   
   procedure Expect (
                     Symbol : in String;
                     Pos : in out Integer;
                     Ch_Set : in Character_Set
                    );
   
   -- Skip the next symbols from the specified character set:
   
   procedure Skip (
                   Symbol : in String;
                   Pos : in out Integer;
                   Ch_Set : in Character_Set
                  );
   
   -- Check whether the string S has only characters from the specified set:
   
   function Has_Only_Characters (S : String; CS : Character_Set) return Boolean;

   -- Check whether the string S has some characters from the specified set:
   
   function Has_Characters (S : String; CS : Character_Set) return Boolean;
   
end Parser_Tools;
