with Text_IO;                   use Text_IO;
with Ada.Integer_Text_IO;       use Ada.Integer_Text_IO;
with Ada.Command_Line;          use Ada.Command_Line;
with Ada.Environment_Variables; use Ada.Environment_Variables;
with Ada.Strings.Fixed;         use Ada.Strings.Fixed;
with Ada.Strings.Maps;          use Ada.Strings.Maps;

with Symmetry_Operations;       use Symmetry_Operations;
with Shmueli_Matrices;          use Shmueli_Matrices;

with Shmueli_Symbol_Parser;     use Shmueli_Symbol_Parser;

with Project_Version;

procedure Decode_Shmueli is
   
   -- This program decodes Shmueli explicit symbols [1,2].
   --
   
   -- [1] Shmueli, U. (1984) Space-group algorithms. I. The space
   --  group and its symmetry elements. Acta Crystallographica Section
   --  A Foundations of Crystallography 40(5), 559-567. International
   --  Union of Crystallography (IUCr). DOI:
   --  https://doi.org/10.1107/s0108767384001161
   --
   -- Table 2. Point-group generators (proper rotations)
   
   -- [2] U. Shmueli (ed.) (2001) International Tables for
   --  Crystallography, Vol. B. Reciprocal Space. (c) International
   --  Union of Crystallography 2001. Published by Kluwer Academic
   --  Publishers, Dordrecht/Boston/London. p. 108.

   
   Debug_Print_Matrices : Boolean := False;
   
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
         Put_Line ("Decode a Shmueli Crystallographic spacegroup symbol " &
                     "given on the command line");
         Put_Line ("and output symmtetry operators as general position " &
                     "point coordinates, e.g. '-X,-Y,Z+1/2'");
         New_Line;
         Put_Line ("USAGE:");
         Put_Line ("  " & Command_Name & " 'CON$P2C006$I2A000'");
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
              Decode_Shmueli_Symbol (Argument (I));
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

end Decode_Shmueli;
