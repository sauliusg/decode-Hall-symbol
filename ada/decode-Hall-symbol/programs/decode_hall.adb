with Text_IO;             use Text_IO;
with Ada.Integer_Text_IO; use Ada.Integer_Text_IO;
with Ada.Command_Line;    use Ada.Command_Line;

procedure Decode_Hall is
   
   type Symop is array (1..4, 1..4) of Float;
   
   type Symop_Array is array (Positive range <>) of Symop;
   
   Zero_Matrix : constant Symop := (others => (others => 0.0));
   
   Unity_Matrix : constant Symop :=
     (
      (1 => 1.0, others => 0.0),
      (2 => 1.0, others => 0.0),
      (3 => 1.0, others => 0.0),
      (4 => 1.0, others => 0.0)
     );
   
   Ci_Matrix : constant Symop :=
     (
      (1 => -1.0, others => 0.0),
      (2 => -1.0, others => 0.0),
      (3 => -1.0, others => 0.0),
      (4 => -1.0, others => 0.0)
     );
   
   function Transpose (S : Symop) return Symop is
      R : Symop;
   begin
      for I in S'Range(1) loop
         for J in S'Range(1) loop
            R (I,J) := S (J,I);
         end loop;
      end loop;
      return R;
   end;
   
   A_Translation : constant Symop :=
     Transpose ((4 => (0.0, 0.5, 0.5, 1.0), others => (others => 0.0)));
   
   B_Translation : constant Symop :=
     Transpose ((4 => (0.5, 0.0, 0.5, 1.0), others => (others => 0.0)));
   
   C_Translation : constant Symop :=
     Transpose ((4 => (0.5, 0.5, 0.0, 1.0), others => (others => 0.0)));
   
   I_Translation : constant Symop :=
     Transpose ((4 => (0.5, 0.5, 0.5, 1.0), others => (others => 0.0)));
   
   R_Translation_1 : constant Symop := 
     Transpose ((
                 4 => (1.0/3.0, 2.0/3.0, 2.0/3.0, 0.0),
                 others => (others => 0.0)
                ));
   
   R_Translation_2 : constant Symop :=
     Transpose ((
                 4 => (2.0/3.0, 1.0/3.0, 1.0/3.0, 0.0), others => (others => 0.0)
                ));
   
   F_Translation_1 : constant Symop := A_Translation;
   F_Translation_2 : constant Symop := B_Translation;
   F_Translation_3 : constant Symop := C_Translation;
   
   Principal_Rotations : constant array (1..3, 1..4) of Symop :=
     (
      -- axis x (a)
      1 => (
            1 => ( -- twofold axis
                  (1.0,  0.0,  0.0,  0.0),
                  (0.0, -1.0,  0.0,  0.0),
                  (0.0,  0.0, -1.0,  0.0),
                  (0.0,  0.0,  0.0,  1.0)
                 ),
            2 => ( --treefold axis
                  (1.0,  0.0,  0.0,  0.0),
                  (0.0,  0.0, -1.0,  0.0),
                  (0.0,  1.0, -1.0,  0.0),
                  (0.0,  0.0,  0.0,  1.0)
                 ),
            3 => ( -- fourfold axis
                  (1.0,  0.0,  0.0,  0.0),
                  (0.0,  0.0, -1.0,  0.0),
                  (0.0,  1.0,  0.0,  0.0),
                  (0.0,  0.0,  0.0,  1.0)
                 ),
            4 => ( --sixfold axis
                  (1.0,  0.0,  0.0,  0.0),
                  (0.0,  1.0, -1.0,  0.0),
                  (0.0,  1.0,  0.0,  0.0),
                  (0.0,  0.0,  0.0,  1.0)
                 )
           ),
      -- axis y (b)
      2 => (
            1 => (
                  (-1.0,  0.0,  0.0,  0.0),
                  ( 0.0,  1.0,  0.0,  0.0),
                  ( 0.0,  0.0, -1.0,  0.0),
                  ( 0.0,  0.0,  0.0,  1.0)
                 ),
            2 => (
                  (-1.0,  0.0,  1.0,  0.0),
                  ( 0.0,  1.0,  0.0,  0.0),
                  (-1.0,  0.0,  0.0,  0.0),
                  ( 0.0,  0.0,  0.0,  1.0)
                 ),
            3 => (
                  ( 0.0,  0.0,  1.0,  0.0),
                  ( 0.0,  1.0,  0.0,  0.0),
                  (-1.0,  0.0,  0.0,  0.0),
                  ( 0.0,  0.0,  0.0,  1.0)
                 ),
            4 => (
                  ( 0.0,  0.0,  1.0,  0.0),
                  ( 0.0,  1.0,  0.0,  0.0),
                  (-1.0,  0.0,  1.0,  0.0),
                  ( 0.0,  0.0,  0.0,  1.0)
                 )
           ),
      -- axis z (c)
      3 => (
            1 => (
                  (-1.0,  0.0,  0.0,  0.0),
                  ( 0.0, -1.0,  0.0,  0.0),
                  ( 0.0,  0.0,  1.0,  0.0),
                  ( 0.0,  0.0,  0.0,  1.0)
                 ),
            2 => (
                  ( 0.0, -1.0,  0.0,  0.0),
                  ( 1.0, -1.0,  0.0,  0.0),
                  ( 0.0,  0.0,  1.0,  0.0),
                  ( 0.0,  0.0,  0.0,  1.0)
                 ),
            3 => (
                  ( 0.0, -1.0,  0.0,  0.0),
                  ( 1.0,  0.0,  0.0,  0.0),
                  ( 0.0,  0.0,  1.0,  0.0),
                  ( 0.0,  0.0,  0.0,  1.0)
                 ),
            4 => (
                  ( 1.0, -1.0,  0.0,  0.0),
                  ( 1.0,  0.0,  0.0,  0.0),
                  ( 0.0,  0.0,  1.0,  0.0),
                  ( 0.0,  0.0,  0.0,  1.0)
                 )
           )

     );
   
   Body_Diagonal_Rotation : constant Symop :=
     (
      (0.0, 0.0, 1.0, 0.0),
      (1.0, 0.0, 0.0, 0.0),
      (0.0, 1.0, 0.0, 0.0),
      (0.0, 0.0, 0.0, 1.0)
     );
   
   procedure Init_Zero (S: out Symop) is
   begin
      S := Zero_Matrix;
   end;
   
   procedure Init_Unity (S: in out Symop) is
   begin
      S := Unity_Matrix;
   end;

   procedure Init_Ci (S: in out Symop) is
   begin
      S := Ci_Matrix;
   end;
   
   procedure Skip_Spaces (S : in String; Pos : in out Integer ) is
   begin
      while Pos <= S'Last and then S (Pos) = ' ' loop
         Pos := Pos + 1;
      end loop;
   end;
   
   procedure Init_Global_Centering
     (
      Symbol : in String;
      Pos : in out Positive;
      Symops : in out Symop_Array;
      N_Symops : in out Positive
     ) is
   begin
      Skip_Spaces (Symbol, Pos);
      if Symbol (Pos) = '-' then
         Init_Ci (Symops (N_Symops));
         N_Symops := N_Symops + 1;
      end if;
   end;
   
   procedure Init_Centering_Matrices 
     (
      Symbol : in String;
      Pos : in out Positive;
      Symops : in out Symop_Array;
      N_Symops : in out Positive
     ) is
   begin
      Skip_Spaces (Symbol, Pos);
      null;
   end;
   
   function Decode_Hall (Symbol : in String) return Symop_Array is
      Max_Symops : constant Integer := 96;
      Symops : Symop_Array (1 .. Max_Symops);
      N_Symops : Positive := 2; -- The first element is allocated for the
                                -- unity matrix.
      Pos : Positive := 1;      -- current position in the string 'Symbol'.
   begin
      Init_Unity (Symops (1));
      Init_Global_Centering (Symbol, Pos, Symops, N_Symops);
      Init_Centering_Matrices (Symbol, Pos, Symops, N_Symops);
      return Symops (1..N_Symops);
   end;
   
   procedure Print_Symop (S : Symop) is
   begin
      for J in Symop'Range(1) loop
         for K in Symop'Range(1) loop
            Put (" " & S (J,K)'Image);
         end loop;
         New_Line;
      end loop;
   end;
   
begin
   
   Print_Symop (A_Translation);
   New_Line;
   Print_Symop (B_Translation);
   New_Line;
   Print_Symop (C_Translation);
   New_Line;
   Print_Symop (I_Translation);
   New_Line;
   New_Line;
   
   Put_Line ("Principal Rotations:");
   for I in Principal_Rotations'Range(1) loop
      for J in Principal_Rotations'Range(2) loop
         Print_Symop (Principal_Rotations (I,J));
         New_Line;         
      end loop;
   end loop;
   
   for I in 1 .. Argument_Count loop
      Put_Line (Argument (I));
      declare
         Symops : Symop_Array := Decode_Hall (Argument (I));
      begin
         for I in Symops'Range loop
            Print_Symop (Symops (I));
            New_Line;
         end loop;
      end;
   end loop;
   
end Decode_Hall;
