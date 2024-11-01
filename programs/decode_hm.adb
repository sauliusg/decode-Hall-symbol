with Text_IO;                   use Text_IO;
with Ada.Integer_Text_IO;       use Ada.Integer_Text_IO;
with Ada.Command_Line;          use Ada.Command_Line;
with Ada.Environment_Variables; use Ada.Environment_Variables;
with Ada.Strings.Fixed;         use Ada.Strings.Fixed;
with Ada.Strings.Maps;          use Ada.Strings.Maps;

with Symmetry_Operations;       use Symmetry_Operations;
with HM_Symbol_Parser;          use HM_Symbol_Parser;

with Project_Version;

procedure Decode_HM is
   
   -- This program decodes Hermann-Mauguin Crystallographic space
   --  group symbols and outputs space group symmetry operators as
   --  lists of general position coordinates.
   
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
              Decode_Hermann_Mauguin_Symbol (Argument (I));
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
