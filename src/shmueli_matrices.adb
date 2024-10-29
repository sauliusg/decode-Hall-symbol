with GCD_Mod; use GCD_Mod;

package body Shmueli_Matrices is
   
   -- "t1 t2 t3 , u1 u2 u3 , v1 v2 v3 – components of the translation
   --  parts of the generators, given in units of 1/12 ; e.g. the
   --  translation part (0 1/2 3/4) is given in Table A1.4.2.1 as 069."
   --  [2, p. 108]
   
   procedure Decode_Shmueli_Symbol_Translation
     (
      Translation_Code : in String;
      Translation      : out Crystallographic_Translation
     ) is
      Code : String (1 .. 3) := Translation_Code;
   begin
      if Translation_Code'Length /= 3 then
         raise CONSTRAINT_ERROR with
           "translation in a Shmuely symbol must be three characters long, " &
           "not " & Integer'Image (Translation_Code'Length) & " as in " &
           "'" & Translation_Code & "'";
      end if;
      for I in Translation'Range loop
         declare
            Value : Natural :=
              Character'Pos (Code (I)) - Character'Pos ('0');
            Divisor : Natural;
         begin
            -- "An exception: (0 0 5/6) is denoted by 005 and not by
            --  0010" [2, p. 108]
            if Value = 5 then
               Value := 10;
            end if;
            if Value /= 0 then
               Divisor := GCD (Value, 12);
               Translation (I) := (Value/Divisor, 12/Divisor);
            else 
               Translation (I) := (0, 1);
            end if;
         end;
      end loop;
   end;
   
end Shmueli_Matrices;
