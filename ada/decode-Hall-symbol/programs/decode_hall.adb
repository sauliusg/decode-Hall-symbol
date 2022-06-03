with Text_IO;             use Text_IO;
with Ada.Integer_Text_IO; use Ada.Integer_Text_IO;
with Ada.Command_Line;    use Ada.Command_Line;

procedure Decode_Hall is
   
   type Symop is array (1..4, 1..4) of Float;
   
   type Symop_Array is array (Positive range <>) of Symop;
   
   procedure Init_Zero (S: out Symop) is
   begin
      for I in S'Range(1) loop
         for J in S'Range(2) loop
            S (I,J) := 0.0;
         end loop;
      end loop;
   end;
   
   procedure Init_Unity (S: in out Symop) is
   begin
      Init_Zero (S);
      for I in S'Range(1) loop
         S (I,I) := 1.0;
      end loop;
   end;

   procedure Init_Ci (S: in out Symop) is
   begin
      Init_Zero (S);
      for I in S'Range(1) loop
         S (I,I) := -1.0;
      end loop;
   end;
   
   procedure Init_Centering_Matrices 
     (
      Symbol : in String;
      Pos : in out Positive;
      Symops : in Symop_Array;
      N_Symops : in out Positive
     ) is
   begin
      null;
   end;
   
   procedure Skip_Spaces (S : in String; Pos : in out Integer ) is
   begin
      while Pos <= S'Last and then S (Pos) = ' ' loop
         Pos := Pos + 1;
      end loop;
   end;
   
   function Decode_Hall (Symbol : in String) return Symop_Array is
      Max_Symops : constant Integer := 96;
      Symops : Symop_Array (1 .. Max_Symops);
      N_Symops : Positive := 1;
      Pos : Positive := 1; -- current position in the string 'Symbol'.
   begin
      Init_Unity (Symops (1));
      Skip_Spaces (Symbol, Pos);
      if Symbol (Pos) = '-' then 
         Init_Ci (Symops (2));
         N_Symops := 2;
      end if;
      Skip_Spaces (Symbol, Pos);
      Init_Centering_Matrices (Symbol, Pos, Symops, N_Symops);
      
      return Symops (1..N_Symops);
   end;
   
begin
   
   for I in 1 .. Argument_Count loop
      Put_Line (Argument (I));
   end loop;
   
end Decode_Hall;
