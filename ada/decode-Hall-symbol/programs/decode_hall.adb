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
   
   type Crystallographic_Translation_Component is record
      Numerator : Integer range 0..6;
      Denominator : Integer range 1..6;
   end record;
      
   type Crystallographic_Translation is array (1..3)
     of Crystallographic_Translation_Component;
   
   A_Translation_Vector : constant Crystallographic_Translation :=
     ((0,1), (1,2), (1,2));
   
   B_Translation_Vector : constant Crystallographic_Translation :=
     ((1,2), (0,1), (1,2));
   
   C_Translation_Vector : constant Crystallographic_Translation :=
     ((1,2), (1,2), (0,1));
   
   I_Translation_Vector : constant Crystallographic_Translation :=
     ((1,2), (1,2), (1,2));
   
   R_Translation_Vector_1 : constant Crystallographic_Translation :=
     ((1,3), (2,3), (2,3));
   
   R_Translation_Vector_2 : constant Crystallographic_Translation :=
     ((2,3), (1,3), (1,3));
   
   F_Translation_Vector_1 : constant Crystallographic_Translation :=
     A_Translation_Vector;
   
   F_Translation_Vector_2 : constant Crystallographic_Translation :=
     B_Translation_Vector;
   
   F_Translation_Vector_3 : constant Crystallographic_Translation :=
     C_Translation_Vector;
   
   
   Translation_a : constant Crystallographic_Translation :=
     ((1,2), (0,1), (0,1));
   Translation_b : constant Crystallographic_Translation :=
     ((0,1), (1,2), (0,1));
   Translation_c : constant Crystallographic_Translation :=
     ((0,1), (0,1), (1,2));
   Translation_n : constant Crystallographic_Translation :=
     ((1,2), (1,2), (1,2));
   Translation_u : constant Crystallographic_Translation :=
     ((1,4), (0,1), (0,1));
   Translation_v : constant Crystallographic_Translation :=
     ((0,1), (1,4), (0,1));
   Translation_w : constant Crystallographic_Translation :=
     ((0,1), (0,1), (1,4));
   Translation_d : constant Crystallographic_Translation :=
     ((1,4), (1,4), (1,4));
   
   Translations_3_1 : constant array (1..3) of Crystallographic_Translation :=
     (((1,3), (0,1), (0,1)), ((0,1), (1,3), (0,1)), ((0,1), (0,1), (1,3)));
   Translations_3_2 : constant array (1..3) of Crystallographic_Translation :=
     (((2,3), (0,1), (0,1)), ((0,1), (2,3), (0,1)), ((0,1), (0,1), (2,3)));
   Translations_4_1 : constant array (1..3) of Crystallographic_Translation :=
     (((1,4), (0,1), (0,1)), ((0,1), (1,4), (0,1)), ((0,1), (0,1), (1,4)));
   Translations_4_3 : constant array (1..3) of Crystallographic_Translation :=
     (((3,4), (0,1), (0,1)), ((0,1), (3,4), (0,1)), ((0,1), (0,1), (3,4)));
   Translations_6_1 : constant array (1..3) of Crystallographic_Translation :=
     (((1,6), (0,1), (0,1)), ((0,1), (1,6), (0,1)), ((0,1), (0,1), (1,6)));
   Translations_6_2 : constant array (1..3) of Crystallographic_Translation :=
     (((2,6), (0,1), (0,1)), ((0,1), (2,6), (0,1)), ((0,1), (0,1), (2,6)));
   Translations_6_4 : constant array (1..3) of Crystallographic_Translation :=
     (((4,6), (0,1), (0,1)), ((0,1), (4,6), (0,1)), ((0,1), (0,1), (4,6)));
   Translations_6_5 : constant array (1..3) of Crystallographic_Translation :=
     (((5,6), (0,1), (0,1)), ((0,1), (5,6), (0,1)), ((0,1), (0,1), (5,6)));
   
   function To_Symop (T : Crystallographic_Translation) return Symop is
      S : Symop := Zero_Matrix;
   begin
      for I in T'Range loop
         S (I,4) := Float (T (I).Numerator) / Float (T (I).Denominator);
      end loop;
      return S;
   end;
   
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
   
   Face_Diagonal_Rotations : constant array (1..3,1..2) of Symop :=
     (
      1 => (
            1 => (
                  (-1.0,  0.0,  0.0,  0.0),
                  ( 0.0,  0.0, -1.0,  0.0),
                  ( 0.0, -1.0,  0.0,  0.0),
                  ( 0.0,  0.0,  0.0,  1.0)
                 ),
            2 => (
                  (-1.0,  0.0,  0.0,  0.0),
                  ( 0.0,  0.0,  1.0,  0.0),
                  ( 0.0,  1.0,  0.0,  0.0),
                  ( 0.0,  0.0,  0.0,  1.0)
                 )
           ),
      2 => (
            1 => (
                  ( 0.0,  0.0, -1.0,  0.0),
                  ( 0.0, -1.0,  0.0,  0.0),
                  (-1.0,  0.0,  0.0,  0.0),
                  ( 0.0,  0.0,  0.0,  1.0)
                 ),
            2 => (
                  ( 0.0,  0.0,  1.0,  0.0),
                  ( 0.0, -1.0,  0.0,  0.0),
                  ( 1.0,  0.0,  0.0,  0.0),
                  ( 0.0,  0.0,  0.0,  1.0)
                 )
           ),
      3 => (
            1 => (
                  ( 0.0, -1.0,  0.0,  0.0),
                  (-1.0,  0.0,  0.0,  0.0),
                  ( 0.0,  0.0, -1.0,  0.0),
                  ( 0.0,  0.0,  0.0,  1.0)
                 ),
            2 => (
                  ( 0.0,  1.0,  0.0,  0.0),
                  ( 1.0,  0.0,  0.0,  0.0),
                  ( 0.0,  0.0, -1.0,  0.0),
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
   
   procedure Put (S : Symop) is
   begin
      for J in Symop'Range(1) loop
         for K in Symop'Range(1) loop
            Put (" " & S (J,K)'Image);
         end loop;
         New_Line;
      end loop;
   end;
   
   procedure Put (S : Crystallographic_Translation) is
      Ratio : Float;
   begin
      for I in S'Range loop
         Ratio := Float (S (I).Numerator) / Float (S (I).Denominator); 
         Put (" " & Ratio'Image);
      end loop;
      New_Line;
   end;
   
begin
   
   Put (A_Translation_Vector);
   New_Line;
   Put (B_Translation_Vector);
   New_Line;
   Put (C_Translation_Vector);
   New_Line;
   Put (I_Translation_Vector);
   New_Line;
   
   Put_Line ("Translation Matrices:");
   for I in Translations_3_1'Range loop
      Put (To_Symop (Translations_3_1 (I)));
      New_Line;
   end loop;
   
   for I in Translations_6_4'Range loop
      Put (To_Symop (Translations_6_4 (I)));
      New_Line;
   end loop;
   
   for I in Translations_6_5'Range loop
      Put (To_Symop (Translations_6_5 (I)));
      New_Line;
   end loop;
   
   Put_Line ("Principal Rotations:");
   for I in Principal_Rotations'Range(1) loop
      for J in Principal_Rotations'Range(2) loop
         Put (Principal_Rotations (I,J));
         New_Line;         
      end loop;
   end loop;
   
   Put_Line ("Face Diagonal Rotations:");
   for I in Face_Diagonal_Rotations'Range(1) loop
      for J in Face_Diagonal_Rotations'Range(2) loop
         Put (Face_Diagonal_Rotations (I,J));
         New_Line;         
      end loop;
   end loop;
   
   for I in 1 .. Argument_Count loop
      Put_Line (Argument (I));
      declare
         Symops : Symop_Array := Decode_Hall (Argument (I));
      begin
         for I in Symops'Range loop
            Put (Symops (I));
            New_Line;
         end loop;
      end;
   end loop;
   
end Decode_Hall;
