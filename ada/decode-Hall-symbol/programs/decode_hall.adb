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
   
   Inversion_Matrices : constant array (1..2) of Symop :=
     (Unity_Matrix, Ci_Matrix);
   
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
   
   Translations_3_1 : constant Crystallographic_Translation_Component := (1,3);
   Translations_3_2 : constant Crystallographic_Translation_Component := (2,3);
   
   Translations_4_1 : constant Crystallographic_Translation_Component := (1,4);
   Translations_4_3 : constant Crystallographic_Translation_Component := (3,4);
   
   Translations_6_1 : constant Crystallographic_Translation_Component := (1,6);
   Translations_6_2 : constant Crystallographic_Translation_Component := (2,6);
   Translations_6_4 : constant Crystallographic_Translation_Component := (4,6);
   Translations_6_5 : constant Crystallographic_Translation_Component := (5,6);
   
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
   
   function To_Symop (T : Crystallographic_Translation) return Symop is
      S : Symop := Zero_Matrix;
   begin
      for I in T'Range loop
         S (I,4) := Float (T (I).Numerator) / Float (T (I).Denominator);
      end loop;
      return S;
   end;
   
   function To_Symop (T : Crystallographic_Translation_Component;
                      Axis : Positive) return Symop is
      S : Symop := Zero_Matrix;
   begin
      S (Axis,4) := Float (T.Numerator) / Float (T.Denominator);
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
   
   procedure Skip_Spaces (S : in String; Pos : in out Integer ) is
   begin
      while Pos <= S'Last and then S (Pos) = ' ' loop
         Pos := Pos + 1;
      end loop;
   end;
   
   procedure Get_Hall_Symbol_Inversions
     (
      Symbol : in String;
      Pos : in out Integer;
      N_Inversions : out Positive
     )
   is
   begin
      Skip_Spaces (Symbol, Pos);
      if Symbol (Pos) = '-' then
         N_Inversions := 2;
         Pos := Pos + 1;
      else
         N_Inversions := 1;
      end if;
   end;
   
   procedure Get_Hall_Symbol_Centerings
     (
      Symbol : in String;
      Pos : in out Positive;
      Centering : out Symop_Array;
      N_Centering : out Positive
     )
   is
   begin
      Skip_Spaces (Symbol, Pos);
      Centering (1) := Zero_Matrix;
      case Symbol (Pos) is
         when 'P' =>
           N_Centering := 1;
         when 'A' =>
           Centering (2) := To_Symop (A_Translation_Vector);
           N_Centering := 2;
         when 'B' =>
           Centering (2) := To_Symop (B_Translation_Vector);
           N_Centering := 2;
         when 'C' =>
           Centering (2) := To_Symop (C_Translation_Vector);
           N_Centering := 2;
         when 'I' =>
           Centering (2) := To_Symop (I_Translation_Vector);
           N_Centering := 2;
         when 'F' =>
           Centering (2) := To_Symop (F_Translation_Vector_1);
           Centering (3) := To_Symop (F_Translation_Vector_2);
           Centering (4) := To_Symop (F_Translation_Vector_3);
           N_Centering := 4;
         when 'R' =>
           Centering (2) := To_Symop (R_Translation_Vector_1);
           Centering (3) := To_Symop (R_Translation_Vector_2);
           N_Centering := 3;
         when others =>
           Centering (1) := Unity_Matrix;
           N_Centering := 1;           
      end case;
      Pos := Pos + 1;
   end;
   
   function Decode_Hall (Symbol : in String) return Symop_Array is
      Max_Symops : constant Integer := 96;
      Symops : Symop_Array (1 .. Max_Symops);
      N_Symops : Positive := 1; -- The first element is allocated for the
                                -- unity matrix.
      Pos : Positive := 1;      -- current position in the string 'Symbol'.
      
      N_Inversions : Positive;
      
      Rotations : Symop_Array (1..3);
      N_Rotations : Positive;
      
      Centering : Symop_Array (1..4);
      N_Centering : Positive;
      
   begin
      Symops (1) := Unity_Matrix;
      
      Get_Hall_Symbol_Inversions (Symbol, Pos, N_Inversions);
      Get_Hall_Symbol_Centerings (Symbol, Pos, Centering, N_Centering);
      
      Put_Line ("Inversions:");
      for I in 1..N_Inversions loop
         Put (Inversion_Matrices (I));
         New_Line;
      end loop;
      
      Put_Line ("Centerings:");
      for I in 1..N_Centering loop
         Put (Centering (I));
         New_Line;
      end loop;
      
      return Symops (1..N_Symops);
   end;
   
begin
   
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
