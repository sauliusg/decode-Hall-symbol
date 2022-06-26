with Text_IO;                   use Text_IO;
with Ada.Integer_Text_IO;       use Ada.Integer_Text_IO;
with Ada.Command_Line;          use Ada.Command_Line;
with Ada.Environment_Variables; use Ada.Environment_Variables;
with Ada.Strings.Fixed;         use Ada.Strings.Fixed;
with Ada.Strings.Maps;          use Ada.Strings.Maps;

procedure Decode_Hall is
   
   -- This program decodes Hall Crystallographic space group symbols
   --  and outputs space group symmetry operators as lists of
   --  general position coordinates.
   --
   -- The description of Hall symbols and decoding tables are taken
   --  from [1]. Extended notation with the origin shift notation is
   --  taken from [2], and the change-of-basis transformations are
   --  described in [3].
   --
   -- Refs.:
   --
   -- 1. Hall, S. R. "Space-group notation with an explicit
   --  origin". Acta Crystallographica Section A, International Union
   --  of Crystallography (IUCr), 1981, 37, 517-525
   --  DOI: https://doi.org/10.1107/s0567739481001228
   --
   -- 2. Sydney R. Hall, Ralf W. Grosse-Kunstleve "Concise Space-Group
   --  Symbols". 1996, URL: https://cci.lbl.gov/sginfo/hall_symbols.html
   --  [accessed: 2022-06-14T15:24+03:00]
   --
   -- 3. International Tables Volume B 1994, Section 1.4. "Symmetry in
   --  reciprocal space". URL: https://onlinelibrary.wiley.com/iucr/itc/B/
   --  [accessed: 2022-06-14T15:35+03:00]
   
   Debug_Print_Matrices : Boolean := False;
   
   UNKNOWN_AXIS : exception;
   UNKNOWN_ROTATION : exception;
   UNKNOWN_CENTERING : exception;
   UNKNOWN_TRANSLATION : exception;
   UNEXPECTED_SYMBOL : exception;
   
   type Axis_Direction_Type is
     (X_AXIS, Y_AXIS, Z_AXIS, UNKNOWN);
   
   subtype Known_Axis_Direction is
     Axis_Direction_Type range X_AXIS .. Z_AXIS;
   
   type Axis_Order_Type is
     (IDENTITY, TWOFOLD, THREEFOLD, FOURFOLD, SIXFOLD, UNKNOWN);
   
   subtype Known_Axis_Order is
     Axis_Order_Type range IDENTITY .. SIXFOLD;
   
   type Symmetry_Operator is array (1..4, 1..4) of Float;
   
   type Symmetry_Operator_Array is
     array (Positive range <>) of Symmetry_Operator;
   
   Zero_Matrix : constant Symmetry_Operator := (others => (others => 0.0));
   
   Unity_Matrix : constant Symmetry_Operator :=
     (
      (1 => 1.0, others => 0.0),
      (2 => 1.0, others => 0.0),
      (3 => 1.0, others => 0.0),
      (4 => 1.0, others => 0.0)
     );
   
   Ci_Matrix : constant Symmetry_Operator :=
     (
      (1 => -1.0, others => 0.0),
      (2 => -1.0, others => 0.0),
      (3 => -1.0, others => 0.0),
      (4 =>  1.0, others => 0.0)
     );
   
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
   
   procedure Put (F : File_Type; S : Symmetry_Operator) is
   begin
      for J in Symmetry_Operator'Range(1) loop
         for K in Symmetry_Operator'Range(1) loop
            Put (F, " " & S (J,K)'Image);
         end loop;
         New_Line (F);
      end loop;
   end;
   
   procedure Put (F : File_Type; S : Crystallographic_Translation) is
      Ratio : Float;
   begin
      for I in S'Range loop
         Ratio := Float (S (I).Numerator) / Float (S (I).Denominator); 
         Put (F, " " & Ratio'Image);
      end loop;
      New_Line (F);
   end;
   
   function To_Symmetry_Operator (T : Crystallographic_Translation) 
                                 return Symmetry_Operator
   is
      S : Symmetry_Operator := Unity_Matrix;
   begin
      for I in T'Range loop
         S (I,4) := Float (T (I).Numerator) / Float (T (I).Denominator);
      end loop;
      return S;
   end;
   
   function Axis_Index (Direction : Known_Axis_Direction) return Positive is
   begin
      case Direction is
         when X_AXIS => return 1;
         when Y_AXIS => return 2;
         when Z_AXIS => return 3;
      end case;
   end Axis_Index;
      
   function To_Symmetry_Operator (T : Crystallographic_Translation_Component;
                                  Axis_Direction : Known_Axis_Direction)
                                 return Symmetry_Operator
   is
      S : Symmetry_Operator := Unity_Matrix;      
   begin
      S (Axis_Index (Axis_Direction), 4) :=
        Float (T.Numerator) / Float (T.Denominator);
      return S;
   end To_Symmetry_Operator;
   
   Principal_Rotations : constant array 
     (Known_Axis_Direction, Known_Axis_Order) of Symmetry_Operator :=
     (
      X_AXIS => 
        ( -- axis x (a)
          IDENTITY =>
            (
             (1.0,  0.0,  0.0,  0.0),
             (0.0,  1.0,  0.0,  0.0),
             (0.0,  0.0,  1.0,  0.0),
             (0.0,  0.0,  0.0,  1.0)                   
            ),
          TWOFOLD =>
            (
             (1.0,  0.0,  0.0,  0.0),
             (0.0, -1.0,  0.0,  0.0),
             (0.0,  0.0, -1.0,  0.0),
             (0.0,  0.0,  0.0,  1.0)
            ),
          THREEFOLD =>
            (
             (1.0,  0.0,  0.0,  0.0),
             (0.0,  0.0, -1.0,  0.0),
             (0.0,  1.0, -1.0,  0.0),
             (0.0,  0.0,  0.0,  1.0)
            ),
          FOURFOLD =>
            (
             (1.0,  0.0,  0.0,  0.0),
             (0.0,  0.0, -1.0,  0.0),
             (0.0,  1.0,  0.0,  0.0),
             (0.0,  0.0,  0.0,  1.0)
            ),
          SIXFOLD =>
            (
             (1.0,  0.0,  0.0,  0.0),
             (0.0,  1.0, -1.0,  0.0),
             (0.0,  1.0,  0.0,  0.0),
             (0.0,  0.0,  0.0,  1.0)
            )
        ),
      Y_AXIS => 
        ( -- axis y (b)
          IDENTITY =>
            (
             (1.0,  0.0,  0.0,  0.0),
             (0.0,  1.0,  0.0,  0.0),
             (0.0,  0.0,  1.0,  0.0),
             (0.0,  0.0,  0.0,  1.0)                   
            ),
          TWOFOLD =>
            (
             (-1.0,  0.0,  0.0,  0.0),
             ( 0.0,  1.0,  0.0,  0.0),
             ( 0.0,  0.0, -1.0,  0.0),
             ( 0.0,  0.0,  0.0,  1.0)
            ),
          THREEFOLD =>
            (
             (-1.0,  0.0,  1.0,  0.0),
             ( 0.0,  1.0,  0.0,  0.0),
             (-1.0,  0.0,  0.0,  0.0),
             ( 0.0,  0.0,  0.0,  1.0)
            ),
          FOURFOLD =>
            (
             ( 0.0,  0.0,  1.0,  0.0),
             ( 0.0,  1.0,  0.0,  0.0),
             (-1.0,  0.0,  0.0,  0.0),
             ( 0.0,  0.0,  0.0,  1.0)
            ),
          SIXFOLD =>
            (
             ( 0.0,  0.0,  1.0,  0.0),
             ( 0.0,  1.0,  0.0,  0.0),
             (-1.0,  0.0,  1.0,  0.0),
             ( 0.0,  0.0,  0.0,  1.0)
            )
        ),
      -- axis z (c)
      Z_AXIS =>
        (
         IDENTITY =>
           ( -- identity
             (1.0,  0.0,  0.0,  0.0),
             (0.0,  1.0,  0.0,  0.0),
             (0.0,  0.0,  1.0,  0.0),
             (0.0,  0.0,  0.0,  1.0)                   
           ),
         TWOFOLD =>
           (
            (-1.0,  0.0,  0.0,  0.0),
            ( 0.0, -1.0,  0.0,  0.0),
            ( 0.0,  0.0,  1.0,  0.0),
            ( 0.0,  0.0,  0.0,  1.0)
           ),
         THREEFOLD =>
           (
            ( 0.0, -1.0,  0.0,  0.0),
            ( 1.0, -1.0,  0.0,  0.0),
            ( 0.0,  0.0,  1.0,  0.0),
            ( 0.0,  0.0,  0.0,  1.0)
           ),
         FOURFOLD =>
           (
            ( 0.0, -1.0,  0.0,  0.0),
            ( 1.0,  0.0,  0.0,  0.0),
            ( 0.0,  0.0,  1.0,  0.0),
            ( 0.0,  0.0,  0.0,  1.0)
           ),
         SIXFOLD
           => 
           (
            ( 1.0, -1.0,  0.0,  0.0),
            ( 1.0,  0.0,  0.0,  0.0),
            ( 0.0,  0.0,  1.0,  0.0),
            ( 0.0,  0.0,  0.0,  1.0)
           )
        )
     );
   
   Face_Diagonal_Rotations :
     constant array (Known_Axis_Direction,1..2) of Symmetry_Operator :=
     (
      X_AXIS => (
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
      Y_AXIS => (
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
      Z_AXIS => (
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
   
   Body_Diagonal_Rotation : constant Symmetry_Operator :=
     (
      (0.0, 0.0, 1.0, 0.0),
      (1.0, 0.0, 0.0, 0.0),
      (0.0, 1.0, 0.0, 0.0),
      (0.0, 0.0, 0.0, 1.0)
     );
   
   function Machine_Epsilon return Float is
      Epsilon : Float := 1.0;
   begin
      while 1.0 + Epsilon /= 1.0 loop
         Epsilon := Epsilon / 2.0;
      end loop;
      return Epsilon;
   end;
   
   Eps : constant Float := 16.0 * Machine_Epsilon;

   procedure Snap_To_Crystallographic_Translations
     (M : in out Symmetry_Operator) is
   begin
      for I in 1 .. M'Last(1) - 1 loop
         M (I,4) := M (I,4) - Float'Floor (M (I,4));
         if abs (M (I,4) - 1.0/3.0) < Eps  then
            M (I,4) := 1.0/3.0;
         elsif abs (M (I,4) - 2.0/3.0) < Eps then
            M (I,4) := 2.0/3.0;
         elsif abs (M (I,4) - 1.0/6.0) < Eps then
            M (I,4) := 1.0/6.0;
         elsif abs (M (I,4) - 5.0/6.0) < Eps then
            M (I,4) := 5.0/6.0;
         elsif abs (M (I,4) - 1.0/2.0) < Eps then
            M (I,4) := 1.0/2.0;
         elsif abs (M (I,4) - 1.0/4.0) < Eps then
            M (I,4) := 1.0/4.0;
         elsif abs (M (I,4) - 3.0/4.0) < Eps then
            M (I,4) := 3.0/4.0;
         elsif abs (M (I,4) - 0.0) < Eps or else abs (M (I,4) - 1.0) < Eps then
            M (I,4) := 0.0;
         end if;
      end loop;
   end Snap_To_Crystallographic_Translations;
   
   procedure Add
     (M : in out Symmetry_Operator; T : Crystallographic_Translation) is
   begin
      for I in 1..3 loop
         M (I,4) := M (I,4) + 
           Float (T (I).Numerator) / Float (T (I).Denominator);
      end loop;
      Snap_To_Crystallographic_Translations (M);
   end;
   
   procedure Add 
     (
      M : in out Symmetry_Operator;
      T : Crystallographic_Translation_Component;
      Axis_Direction : Known_Axis_Direction
     ) 
   is
      I : Positive := Axis_Index (Axis_Direction);
   begin
      M(I,4) := M(I,4) +
        Float (T.Numerator) / Float (T.Denominator);
      Snap_To_Crystallographic_Translations (M);
   end;
   
   function "*" (M1, M2 : Symmetry_Operator) return Symmetry_Operator is
      M : Symmetry_Operator;
   begin
      pragma Assert (M1'Last(2) = M2'Last(1));
      
      for I in M1'Range(1) loop
         for J in M2'Range(2) loop
            M (I,J) := 0.0;
            for K in M2'Range(1) loop
               M (I,J) := M (I,J) + M1 (I,K) * M2(K,J);
            end loop;
         end loop;
      end loop;
      Snap_To_Crystallographic_Translations (M);
      return M;
   end "*";
   
   type Matrix3x3 is array (1..3,1..3) of Float;
   
   function Det (M : Matrix3x3) return Float is
   ( 
     M(1,1)*M(2,2)*M(3,3) + 
       M(1,2)*M(2,3)*M(3,1) + 
       M(2,1)*M(3,2)*M(1,3) -
       M(1,3)*M(2,2)*M(3,1) - 
       M(2,1)*M(1,2)*M(3,3) - 
       M(1,1)*M(2,3)*M(3,2)
   );
      
   function Invert (M : Matrix3x3) return Matrix3x3 is
      D : Float := Det (M);
      
      function Adjunct (P, Q : Integer) return Float is
         A : array (1..4,1..2) of Float;
         K, L : Integer;
         Coeff : Float;
      begin
         K := 1; L := 1;
         for I in M'Range(1) loop
            if I /= P then 
               for J in M'Range(2) loop
                  if J /= Q then
                     A(K,L) := M(I,J);
                     L := L + 1;
                  end if;
               end loop;
               L := 1;
               K := K + 1;
            end if;
         end loop;
         if (P + Q) mod 2 = 0 then
            Coeff := 1.0;
         else
            Coeff := -1.0;
         end if;
         return Coeff * (A(1,1)*A(2,2) - A(1,2)*A(2,1));
      end Adjunct;
      
   begin -- Invert
      return (
              ( Adjunct(1,1)/D, Adjunct(2,1)/D, Adjunct(3,1)/D ),
              ( Adjunct(1,2)/D, Adjunct(2,2)/D, Adjunct(3,2)/D ),
              ( Adjunct(1,3)/D, Adjunct(2,3)/D, Adjunct(3,3)/D )
             );
   end Invert;
   
   -- A very specific Determinant routine for symops:
   function Det (S : Symmetry_Operator) return Float is
      R : Matrix3x3;
   begin
      for I in R'Range(1) loop
         for J in R'Range(2) loop
            R (I,J) := S (I,J);
         end loop;
      end loop;
      return Det (R);
   end;
   
   -- a very specific inversion routine for sympos:
   function Invert (S : Symmetry_Operator) return Symmetry_Operator is
      R : Matrix3x3;
      Inv : Symmetry_Operator;
   begin
      for I in R'Range(1) loop
         for J in R'Range(2) loop
            R (I,J) := S (I,J);
         end loop;
      end loop;
      
      -- rotation matrix of the inverse symop is an inverted rotation
      --  matrix:
      R := Invert (R);
      
      for I in R'Range(1) loop
         for J in R'Range(2) loop
            Inv (I,J) := R (I,J);
         end loop;
      end loop;
      
      -- compute the inverse translation:
      
      -- assume R' is the inverse of R:
      -- R * R' = I, where I is the unity matrix.
      -- then:
      -- (R,t) * (R',t') = (R*R', R't + t') = (I,0)
      -- =>
      -- R't + t' = 0
      -- =>
      -- t' = -R' * t
      
      for I in R'Range(1) loop
         Inv (I,4) := 0.0;
         for J in R'Range(2) loop
            Inv (I,4) := Inv (I,4) - R (I,J) * S (J,4);
         end loop;
      end loop;
      
      -- The last row is (0,0,0,1):
      Inv (4,1) := 0.0;
      Inv (4,2) := 0.0;
      Inv (4,3) := 0.0;
      Inv (4,4) := 1.0;
      
      return Inv;
   end Invert;
   
   -- -------------------------------------------------------------------------
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
   end Get_Hall_Symbol_Inversions;
   
   procedure Get_Hall_Symbol_Centerings
     (
      Symbol : in String;
      Pos : in out Positive;
      Centering : out Symmetry_Operator_Array;
      N_Centering : out Positive
     )
   is
   begin
      Skip_Spaces (Symbol, Pos);
      Centering (1) := Unity_Matrix;
      case Symbol (Pos) is
         when 'P' =>
           N_Centering := 1;
         when 'A' =>
           Centering (2) := To_Symmetry_Operator (A_Translation_Vector);
           N_Centering := 2;
         when 'B' =>
           Centering (2) := To_Symmetry_Operator (B_Translation_Vector);
           N_Centering := 2;
         when 'C' =>
           Centering (2) := To_Symmetry_Operator (C_Translation_Vector);
           N_Centering := 2;
         when 'I' =>
           Centering (2) := To_Symmetry_Operator (I_Translation_Vector);
           N_Centering := 2;
         when 'F' =>
           Centering (2) := To_Symmetry_Operator (F_Translation_Vector_1);
           Centering (3) := To_Symmetry_Operator (F_Translation_Vector_2);
           Centering (4) := To_Symmetry_Operator (F_Translation_Vector_3);
           N_Centering := 4;
         when 'R' =>
           Centering (2) := To_Symmetry_Operator (R_Translation_Vector_1);
           Centering (3) := To_Symmetry_Operator (R_Translation_Vector_2);
           N_Centering := 3;
         when others =>
            raise UNKNOWN_CENTERING with
              "unknown centering symbol " & Symbol (Pos)'Image;
      end case;
      Pos := Pos + 1;
   end Get_Hall_Symbol_Centerings;
   
   function Rotation_Axis_Index (Rotation_Character : Character)
                                return Known_Axis_Order
   is
   begin
      case Rotation_Character is
         when '1' => return IDENTITY;
         when '2' => return TWOFOLD;
         when '3' => return THREEFOLD;
         when '4' => return FOURFOLD;
         when '6' => return SIXFOLD;
         when others => 
            raise UNKNOWN_ROTATION
              with "rotation " & Rotation_Character'Image;
      end case;
   end Rotation_Axis_Index;
   
   procedure Get_Rotation_Matrix_From_Axis_And_Rotation
     (
      Matrix : out Symmetry_Operator;
      Axis : in Character;
      Rotation : in Character;
      Preceeding_Axis_Direction : in out Axis_Direction_Type;
      Preceeding_Axis_Order : in out Axis_Order_Type;
      Axis_Direction : out Axis_Direction_Type;
      Axis_Number : in Natural
     )
   is
      Current_Axis_Order : constant Known_Axis_Order := 
        Rotation_Axis_Index (Rotation);
   begin
      Axis_Direction := (if Axis_Number = 1 then Z_AXIS else X_AXIS);
      case Axis is 
         when 'x' =>
            Matrix :=
              Principal_Rotations (X_AXIS, Current_Axis_Order);
            Axis_Direction := X_AXIS;
            Preceeding_Axis_Direction := Axis_Direction;
         when 'y' =>
            Matrix :=
              Principal_Rotations (Y_AXIS, Current_Axis_Order);
            Axis_Direction := Y_AXIS;
            Preceeding_Axis_Direction := Axis_Direction;
         when 'z' =>
            Matrix :=
              Principal_Rotations (Z_AXIS, Current_Axis_Order);
            Axis_Direction := Z_AXIS;
            Preceeding_Axis_Direction := Axis_Direction;
         when ''' =>
            Matrix :=
              Face_Diagonal_Rotations (Preceeding_Axis_Direction, 1);
            Preceeding_Axis_Direction := Axis_Direction;
         when '"' =>
            Matrix :=
              Face_Diagonal_Rotations (Preceeding_Axis_Direction, 2);
            Preceeding_Axis_Direction := Axis_Direction;
         when '*' =>
            Matrix :=
              Body_Diagonal_Rotation;
            Preceeding_Axis_Direction := 
              Known_Axis_Direction'Val (Axis_Number - 1);
         when others =>
            raise UNKNOWN_AXIS with "axis character " & Axis'Image;
      end case;
      
      Preceeding_Axis_Order := Current_Axis_Order;
   end Get_Rotation_Matrix_From_Axis_And_Rotation;
   
   procedure Add_Translation_To_The_Rotation_Matrix
     (
      Matrix : out Symmetry_Operator;
      Rotation : Character;
      Translations : String;
      Axis_Direction : in Axis_Direction_Type
     )
   is
   begin
      for Tr of Translations loop
         case Tr is
            when ' ' => null;
            when 'a' => Add (Matrix, Translation_a);
            when 'b' => Add (Matrix, Translation_b);
            when 'c' => Add (Matrix, Translation_c);
            when 'd' => Add (Matrix, Translation_d);
            when 'u' => Add (Matrix, Translation_u);
            when 'v' => Add (Matrix, Translation_v);
            when 'w' => Add (Matrix, Translation_w);
            when 'n' => Add (Matrix, Translation_n);
            when '1' => 
               case Rotation is 
                  when '3' => 
                     Add (Matrix, Translations_3_1, Axis_Direction);
                  when '4' => 
                     Add (Matrix, Translations_4_1, Axis_Direction);
                  when '6' => 
                     Add (Matrix, Translations_6_1, Axis_Direction);
                  when others =>
                     raise UNKNOWN_ROTATION
                       with "mismatching translation " & 
                       Tr'Image & " for rotation " & Rotation'Image;
               end case;
            when '2' => 
               case Rotation is 
                  when '3' => 
                     Add (Matrix, Translations_3_2, Axis_Direction);
                  when '6' => 
                     Add (Matrix, Translations_6_2, Axis_Direction);
                  when others =>
                     raise UNKNOWN_ROTATION
                       with "mismatching translation " & 
                       Tr'Image & " for rotation " & Rotation'Image;
               end case;
            when '3' => 
               case Rotation is 
                  when '4' => 
                     Add (Matrix, Translations_4_3, Axis_Direction);
                  when others =>
                     raise UNKNOWN_ROTATION
                       with "mismatching translation " & 
                       Tr'Image & " for rotation " & Rotation'Image;
               end case;
            when '4' => 
               case Rotation is 
                  when '6' => 
                     Add (Matrix, Translations_6_4, Axis_Direction);
                  when others =>
                     raise UNKNOWN_ROTATION
                       with "mismatching translation " & 
                       Tr'Image & " for rotation " & Rotation'Image;
               end case;
            when '5' => 
               case Rotation is 
                  when '6' => 
                     Add (Matrix, Translations_6_5, Axis_Direction);
                  when others =>
                     raise UNKNOWN_ROTATION
                       with "mismatching translation " & 
                       Tr'Image & " for rotation " & Rotation'Image;
               end case;
            when others =>
               raise UNKNOWN_TRANSLATION
                 with "translation character " & Tr'Image;
         end case;
      end loop;
   end Add_Translation_To_The_Rotation_Matrix;
   
   procedure Construct_Rotation_Matrix
     (
      Matrix : out Symmetry_Operator;
      Inversion : in Character;
      Axis : in Character;
      Rotation : in Character;
      Translations : String;
      Preceeding_Axis_Direction : in out Axis_Direction_Type;
      Preceeding_Axis_Order : in out Axis_Order_Type;
      Axis_Number : in Natural
     )
   is
      Axis_Direction : Axis_Direction_Type;
   begin      
      Get_Rotation_Matrix_From_Axis_And_Rotation
        ( 
          Matrix, Axis, Rotation,
          Preceeding_Axis_Direction,
          Preceeding_Axis_Order,
          Axis_Direction,
          Axis_Number
        );
      
      if Inversion = '-' then
         Matrix := Ci_Matrix * Matrix;
      end if;
      
      Add_Translation_To_The_Rotation_Matrix
        (
         Matrix, Rotation,
         Translations,
         Axis_Direction
        );
   end Construct_Rotation_Matrix;
      
   function Has_Symmetry_Operator
     (
      Symmetry_Operators : Symmetry_Operator_Array;
      Last_Symmetry_Operator_Index : Positive;
      Lookup_Symmetry_Operator : Symmetry_Operator
     )
     return Boolean is
   begin
      for I in 1 .. Last_Symmetry_Operator_Index loop
         if Symmetry_Operators (I) = Lookup_Symmetry_Operator then
            return True;
         end if;
      end loop;
      return False;
   end;
      
   procedure Get_Inversion_Character
     (
      Symbol : in String;
      Pos : in out Positive;
      Inversion : out Character
     )
   is
   begin
      if Pos <= Symbol'Last and then
        Symbol (Pos) = '-' then
         Inversion := Symbol (Pos);
         Pos := Pos + 1;
      else
         Inversion := ' ';
      end if;
   end;
   
   procedure Get_Rotation_Character
     (
      Symbol : in String;
      Pos : in out Positive;
      Rotation : out Character
     )
   is
   begin
      if Pos <= Symbol'Last then
         pragma Assert (Symbol (Pos) in '2'..'6');
         Rotation := Symbol (Pos);
         Pos := Pos + 1;
      else
         Rotation := ' ';
      end if;
   end;
   
   procedure Get_Axis_Character
     (
      Symbol : in String;
      Pos : in out Positive;
      Axis : out Character;
      Axis_Number : in Positive;
      Rotation_Character : in Character;
      Preceeding_Axis_Order : in Axis_Order_Type
     )
   is
   begin
      if Pos <= Symbol'Last and then
        (
         Symbol (Pos) in 'x'..'z' or else
           Symbol (Pos) = ''' or else
           Symbol (Pos) = '*' or else
           Symbol (Pos) = '"'
        )
      then
         Axis := Symbol (Pos);
         Pos := Pos + 1;
      else
         case Axis_Number is
            when 1 => Axis := 'z';
            when 2..3 =>
               case Rotation_Character is
                  when ' ' => null;
                  when '1' => 
                     Axis := 'z';
                  when '2' => 
                     if Preceeding_Axis_Order = TWOFOLD or else
                       Preceeding_Axis_Order = FOURFOLD then
                        Axis := 'x';
                     elsif Preceeding_Axis_Order = THREEFOLD or else 
                       Preceeding_Axis_Order = SIXFOLD then
                        Axis := ''';
                     else
                        raise UNKNOWN_AXIS with
                          "can not determine rotation axis for " &
                          "the preceeding axis " & 
                          Preceeding_Axis_Order'Image &
                          " and current axis " & 
                          Rotation_Character'Image;
                     end if;
                  when '3' => 
                     Axis := '*';
                  when others =>
                     raise UNKNOWN_ROTATION 
                       with "wrong rotation character " & 
                       Rotation_Character'Image &
                       " for axis number" & Axis_Number'Image;
               end case;
            when 4 =>
               Axis := 'x';
            when others =>
               raise UNKNOWN_AXIS with "axis number" & Axis_Number'Image;
         end case;
      end if;
   end Get_Axis_Character;
   
   procedure Get_Translation_Characters
     (
      Symbol : in String;
      Pos : in out Positive;
      Translations : in out String;
      N_Translations: in out Natural
     ) 
   is
      I : Natural := N_Translations;
   begin
      while Pos <= Symbol'Last and then
        (
         Symbol (Pos) in 'a' .. 'd' or else
           Symbol (Pos) in 'u' .. 'w' or else
           Symbol (Pos) in '1' .. '5' or else
           Symbol (Pos) = 'n' 
        ) loop
         I := I + 1;
         Translations (I) := Symbol (Pos);
         Pos := Pos + 1;
      end loop;
      N_Translations := I;
   end Get_Translation_Characters;
      
   procedure Get_Hall_Symbol_Rotation
     (
      Symbol : in String;
      Pos : in out Positive;
      Rotations : out Symmetry_Operator_Array;
      N_Rotations : in out Natural;
      Preceeding_Axis_Direction : in out Axis_Direction_Type;
      Preceeding_Axis_Order : in out Axis_Order_Type;
      Axis_Number : in Positive
     )
   is      
      Inversion : Character;
      Rotation : Character;
      Axis : Character;
      Translations : String (1..3) := (others => ' ');
      N_Translations : Natural := 0;      
   begin
      Skip_Spaces (Symbol, Pos);
      if Pos <= Symbol'Length and then Symbol (Pos) /= '(' then
         Get_Inversion_Character (Symbol, Pos, Inversion);
         Get_Rotation_Character (Symbol, Pos, Rotation);
         
         Get_Translation_Characters (Symbol, Pos, Translations, N_Translations);
         
         Get_Axis_Character (Symbol, Pos, Axis, Axis_Number, Rotation,
                             Preceeding_Axis_Order);
         
         Get_Translation_Characters (Symbol, Pos, Translations, N_Translations);
         
         if Rotation /= ' ' and then Axis /= ' ' then
            N_Rotations := N_Rotations + 1;
            Construct_Rotation_Matrix (Rotations (N_Rotations), 
                                       Inversion, Axis, Rotation,
                                       Translations, Preceeding_Axis_Direction,
                                       Preceeding_Axis_Order,
                                       Axis_Number);
            
            if Has_Symmetry_Operator( Rotations, N_Rotations - 1, 
                                      Rotations (N_Rotations)) then
               N_Rotations := N_Rotations - 1;
            end if;
         end if;
      end if;
   end Get_Hall_Symbol_Rotation;
   
   subtype Character_Set is Ada.Strings.Maps.Character_Set;
   
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
              "is not expected at position" & Pos'Image &
              " in """ & Symbol & """" &
              ", expecing one of """ &
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

   ----------------------------------------------------------------------------
   -- A simple recursive descent parser for the change-of-basis operator:
   
   -- parse a fractional number (e.g. "2/3") and return its value as Float:
   procedure Inc (Result : in out Float; D : in Float) is
   begin
      Result := Result + D;
   end;
   
   function Get_Number (Symbol : in String; Pos : in out Integer) return Float
   is
      Numerator : Natural := 1;
      Denominator : Natural := 1;
      Final_Pos : Integer := Pos;
   begin
      while Final_Pos <= Symbol'Last and then 
        Symbol (Final_Pos) in '0'..'9'
      loop
         Final_Pos := Final_Pos + 1;
      end loop;
      Numerator := Integer'Value (Symbol (Pos..Final_Pos-1));
      Skip_Spaces (Symbol, Final_Pos);
      Pos := Final_Pos;
      if Pos <= Symbol'Last and then Symbol (Pos) = '/' then
         Pos := Pos + 1;
         Skip_Spaces (Symbol, Pos);
         Final_Pos := Pos;
         while Final_Pos <= Symbol'Last and then 
           Symbol (Final_Pos) in '0'..'9'
         loop
            Final_Pos := Final_Pos + 1;
         end loop;
         Denominator := Integer'Value (Symbol (Pos..Final_Pos-1));
         Pos := Final_Pos;
      end if;
      return Float (Numerator) / Float (Denominator);
   end Get_Number;
   
   -- parse the '+1/2*x' factor:
   procedure Parse_Factor
     (
      Symbol : in String;
      Pos : in out Integer;
      Change_Of_Basis : out Symmetry_Operator;
      Row : in Integer;
      Factor : in Float
     ) is
   begin
      Skip_Spaces (Symbol, Pos);
      Expect (Symbol, Pos, To_Set ("xXyYzY"));
      case Symbol (Pos) is
         when 'x'|'X' => Change_Of_Basis (Row, 1) := Factor;
         when 'y'|'Y' => Change_Of_Basis (Row, 2) := Factor;
         when 'z'|'Z' => Change_Of_Basis (Row, 3) := Factor;
         when others =>
            raise UNEXPECTED_SYMBOL with
              "unexpected character " & Symbol (Pos)'Image;
      end case;
      Pos := Pos + 1;
   end Parse_Factor;
   
   -- Parse the "+x", "y", "-z", "1/2" parts in the "+x-y*1/2:
   procedure Parse_Term
     (
      Symbol : in String;
      Pos : in out Integer;
      Change_Of_Basis : out Symmetry_Operator;
      Row : in Integer
     ) is
      Factor : Float := 1.0;
   begin
      Skip_Spaces (Symbol, Pos);
      
      if Pos <= Symbol'Last then
         if Symbol (Pos) = '+' then
            Pos := Pos + 1;
         elsif Symbol (Pos) = '-' then
            Factor := -1.0;
            Pos := Pos + 1;
         end if;
      end if;
      
      Expect (Symbol, Pos, To_Set ("0123456789xXyYzZ"));
      
      if Pos <= Symbol'Last then
         case Symbol (Pos) is
            when 'x'|'X' => 
               Change_Of_Basis (Row, 1) := Factor;
               Pos := Pos + 1;
            when 'y'|'Y' => 
               Change_Of_Basis (Row, 2) := Factor;
               Pos := Pos + 1;
            when 'z'|'Z' => 
               Change_Of_Basis (Row, 3) := Factor;
               Pos := Pos + 1;
            when '0'..'9' =>
               Factor := Factor * Get_Number (Symbol, Pos);
               Skip_Spaces (Symbol, Pos);
               if Pos <= Symbol'Length and then Symbol (Pos) = '*' then
                  Pos := Pos + 1;
                  Parse_Factor (Symbol, Pos, Change_Of_Basis, Row, Factor);
               else
                  Inc (Change_Of_Basis (Row, 4), Factor);
               end if;
            when others =>
               raise UNEXPECTED_SYMBOL with
                 "unexpected symbol " & Character'Image (Symbol (Pos)) &
                 " in the symop """ & Symbol & """";
         end case;
      end if;
   end Parse_Term;
   
   -- parse the "-x+y*1/2" part in the "-x+y*1/2,-z,y+2/3" operator:
   procedure Parse_Symmetry_Operator_Component
     (
      Symbol : in String;
      Pos : in out Integer;
      Change_Of_Basis : out Symmetry_Operator;
      Row : Integer
     ) is
   begin
      loop
         Parse_Term (Symbol, Pos, Change_Of_Basis, Row);
         Skip_Spaces (Symbol, Pos);
         if Pos > Symbol'Length or else 
           Is_In (Symbol (Pos), To_Set (",)")) then
            exit;
         end if;
      end loop;
   end Parse_Symmetry_Operator_Component;
   
   procedure Interpret_Change_Of_Basis_Matrix
      (
       Symbol : in String;
       Pos : in out Integer;
       Change_Of_Basis : out Symmetry_Operator
      )
   is
   begin
      Change_Of_Basis := Zero_Matrix;
      Change_Of_Basis (4,4) := 1.0;
      Skip (Symbol, Pos, To_Set('('));
      Parse_Symmetry_Operator_Component (Symbol, Pos, Change_Of_Basis, 1);
      Skip (Symbol, Pos, To_Set(','));
      Parse_Symmetry_Operator_Component (Symbol, Pos, Change_Of_Basis, 2);
      Skip (Symbol, Pos, To_Set(','));
      Parse_Symmetry_Operator_Component (Symbol, Pos, Change_Of_Basis, 3);
      Skip (Symbol, Pos, To_Set(')'));
   end Interpret_Change_Of_Basis_Matrix;
   
   procedure Get_Shift_Of_Origin (
                                  Symbol : in String;
                                  Pos : in out Integer;
                                  Change_Of_Basis : out Symmetry_Operator
                                 )
   is
      Shift : Integer;
      Sign : Integer := 1;
      S : Symmetry_Operator := Unity_Matrix;
   begin
      Skip_Spaces (Symbol, Pos);
      
      if Pos <= Symbol'Last and then Symbol (Pos) = '(' then
         Pos := Pos + 1;
         
         for I in 1 .. 3 loop
            Expect (Symbol, Pos, To_Set ("-0123456789"));
            if Symbol (Pos) = '-' then
               Sign := -1;
               Pos := Pos + 1;
            end if;
            Get (Symbol (Pos..Symbol'Last), Shift, Pos);
            Pos := Pos + 1;
            S (I,4) := Float (Sign * Shift) / 12.0;
         end loop;
         
         Expect (Symbol, Pos, To_Set (')'));         
      end if;
      Change_Of_Basis := S;
   end Get_Shift_Of_Origin;
   
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
   
   procedure Get_Change_Of_Basis (
                                  Symbol : in String;
                                  Pos : in out Integer;
                                  Change_Of_Basis : out Symmetry_Operator
                                 )
   is
   begin
      if Has_Only_Characters (Symbol (Pos..Symbol'Last),
                              To_Set ("-( 0123456789)")) then
         Get_Shift_Of_Origin (Symbol, Pos, Change_Of_Basis);
      else
         Interpret_Change_Of_Basis_Matrix (Symbol, Pos, Change_Of_Basis);
      end if;
   end;
   
   procedure Build_Group
     (
      Operators : in out Symmetry_Operator_Array;
      N_Operators : in out Positive
     )
   is
      N, M : Positive := N_Operators;
      New_Operator : Symmetry_Operator;
   begin
      loop
         for I in 2 .. N loop
            New_Operator := Operators (I) * Operators (N);
            if not Has_Symmetry_Operator (Operators, M, New_Operator) then
               M := M + 1;
               Operators (M) := New_Operator;
            end if;
         end loop;
         N := N + 1;
         exit when N > M or else M >= Operators'Last;
      end loop;
      N_Operators := M;
   end Build_Group;
   
   function Decode_Hall (Symbol : in String) return Symmetry_Operator_Array is
      Max_Symmetry_Operators : constant Integer := 192;
      
      Symmetry_Operators :
        Symmetry_Operator_Array (1 .. Max_Symmetry_Operators);
      N_Symmetry_Operators : Positive := 1;
      
      Pos : Positive := 1; -- current position in the string 'Symbol'.
      
      Inversions : array (1..2) of Symmetry_Operator :=
        (Unity_Matrix, Ci_Matrix);
      N_Inversions : Positive;
      
      Max_Centering : constant Positive := 8;
      Centering : Symmetry_Operator_Array (1..Max_Centering);
      N_Centering : Positive;
      
      Preceeding_Axis_Direction : Axis_Direction_Type := UNKNOWN;
      Preceeding_Axis_Order : Axis_Order_Type := UNKNOWN;
      
      Change_Of_Basis : Symmetry_Operator;
      
   begin
      Symmetry_Operators (1) := Unity_Matrix;
      
      Get_Hall_Symbol_Inversions (Symbol, Pos, N_Inversions);
      Get_Hall_Symbol_Centerings (Symbol, Pos, Centering, N_Centering);
      
      for Axis_Number in 1..4 loop
         Get_Hall_Symbol_Rotation (Symbol, Pos, Symmetry_Operators,
                                   N_Symmetry_Operators,
                                   Preceeding_Axis_Direction,
                                   Preceeding_Axis_Order, Axis_Number);
      end loop;
      
      Get_Change_Of_Basis (Symbol, Pos, Change_Of_Basis);
      
      -- Apply the change-of-basis operator:
      
      if Change_Of_Basis /= Unity_Matrix then
         declare
            V : Symmetry_Operator := Change_Of_Basis;
            V_Inv : Symmetry_Operator;
         begin
            V_Inv := Invert (V);
            
            for I in 2..N_Symmetry_Operators loop
               Symmetry_Operators (I) := V * Symmetry_Operators (I) * V_Inv;
            end loop;
            
            for I in 2..N_Centering loop
               Centering (I) := V * Centering (I) * V_Inv;
            end loop;
            
            while N_Centering > 1 and then 
              Centering (N_Centering) = Centering (1) 
            loop
               -- remove centerings that became unit matrices after C-o-B:
               N_Centering := N_Centering - 1;
            end loop;
            
            if N_Inversions = 2 then
               Inversions (2) := V * Inversions (2) * V_Inv;
            end if;
         end;
      end if;
      
      -- Generate additional centering operators:
      
      declare
         type Vector_Components is array (1..4) of Float;
         type Vector_Type is record
            Value : Vector_Components;
         end record;
         
         function "*" (S : Symmetry_Operator; T : Vector_Type)
                      return Vector_Type 
         is
            R : Vector_Type := (Value => (others => 0.0));
         begin
            for I in R.Value'Range loop
               for K in T.Value'Range loop
                  R.Value (I) := R.Value (I) +
                    S (I,K) * T.Value (K);
               end loop;
            end loop;
            return R;
         end;
         
         function To_Symmetry_Operator (T : Vector_Type)
                                       return Symmetry_Operator
         is
            S : Symmetry_Operator := Unity_Matrix;
         begin
            for I in 1..3 loop
               S (I,4) := T.Value (I);
            end loop;
            Snap_To_Crystallographic_Translations (S);
            return S;
         end;
         
         function Is_Centering (T : Vector_Type) return Boolean is            
            function Fract (X : Float) return Float is (X - Float'Floor (X));
         begin
            for Component of T.Value loop
               if abs (Fract (Component)) >= Eps then
                  return True;
               end if;
            end loop;
            return False;
         end;
         
         Unit_Vectors : array (1..3) of Vector_Type :=
           (
            (Value => (1.0, 0.0, 0.0, 1.0)),
            (Value => (0.0, 1.0, 0.0, 1.0)),
            (Value => (0.0, 0.0, 1.0, 1.0))
           );
         
         C_O_B_Rotation : Symmetry_Operator := Change_Of_Basis;
         
      begin
         for I in 1..3 loop
            C_O_B_Rotation (I,4) := 0.0;
         end loop;
         for Vector of Unit_Vectors loop
            if Is_Centering (C_O_B_Rotation * Vector) then
               N_Centering := N_Centering + 1;
               Centering (N_Centering) :=
                 To_Symmetry_Operator (C_O_B_Rotation * Vector);
            end if;
         end loop;
         
         -- see if multiplication of two new centerings gives a third one:
         Build_Group (Centering, N_Centering);
      end;
         
      -- Print out all matrices if requested:
      
      if Debug_Print_Matrices then
         Put_Line (Standard_Error, "Inversions:");
         for I in 1..N_Inversions loop
            Put (Standard_Error, Inversions (I));
            New_Line (Standard_Error);
         end loop;
         
         Put_Line (Standard_Error, "Centerings:");
         for I in 1..N_Centering loop
            Put (Standard_Error, Centering (I));
            New_Line (Standard_Error);
         end loop;
         
         Put_Line (Standard_Error, "Rotations:");
         for I in 1..N_Symmetry_Operators loop
            Put_Line (Standard_Error, I'Image);
            Put (Standard_Error, Symmetry_Operators (I));
            New_Line (Standard_Error);
         end loop;
         
         Put_Line (Standard_Error, "Change of basis:");
         Put (Standard_Error, Change_Of_Basis);
         New_Line (Standard_Error);
      end if;
      
      -- Reconstruct all rotation operators:
      
      Build_Group (Symmetry_Operators, N_Symmetry_Operators);
      
      -- Add centering and inversion matrices:
      
      declare
         M : Positive := N_Symmetry_Operators;
         New_Symmetry_Operator : Symmetry_Operator;
      begin
         for I in 1..N_Inversions loop
            for C in 1..N_Centering loop
               if I /= 1 or else C /= 1 then
                  for S in 1..N_Symmetry_Operators loop
                     New_Symmetry_Operator :=
                       Symmetry_Operators (S) * Centering (C) * Inversions (I);
                     if not Has_Symmetry_Operator (Symmetry_Operators, M, 
                                                   New_Symmetry_Operator) then
                       M := M + 1;
                       Symmetry_Operators (M) := New_Symmetry_Operator;
                     end if;
                  end loop;
               end if;
            end loop;
         end loop;
         N_Symmetry_Operators := M;
      end;      
      
      return Symmetry_Operators (1..N_Symmetry_Operators);
   end Decode_Hall;
   
   function As_String (S : Symmetry_Operator) return String is
      Buffer : String (1..100); -- large enough to hold any symop
      Pos : Positive := 1;
      Non_Zero_Printed : Boolean;
      
      function Rational_Translation (T : Float) return String is
         Buffer : String := T'Image;
         Idx : Positive := 1;
      begin
         if T = 0.0 then
            return "";
         elsif T = 0.5 then
            return "+1/2";
         elsif T = 1.5 then
            return "+3/2";
         elsif T = 1.0/3.0 then
            return "+1/3";
         elsif T = 2.0/3.0 then
            return "+2/3";
         elsif T = 1.0/4.0 then
            return "+1/4";
         elsif T = 3.0/4.0 then
            return "+3/4";
         elsif T = 1.0/6.0 then
            return "+1/6";
         elsif T = 5.0/6.0 then
            return "+5/6";
         elsif T = 1.0/8.0 then
            return "+1/8";
         elsif T = 3.0/8.0 then
            return "+3/8";
         elsif T = 5.0/8.0 then
            return "+5/8";
         elsif T = 7.0/8.0 then
            return "+7/8";
         else
            while Idx <= Buffer'Last and then
              (Buffer (Idx) = ' ' or else
                 Buffer (Idx) = '+' or else
                 Buffer (Idx) = '-') loop
               Idx := Idx + 1;
            end loop;
            return (if T < 0.0 then "-" else "+") & Buffer (Idx..Buffer'Last);
         end if;
      end Rational_Translation;
      
   begin
      for I in 1 .. S'Last(2) - 1 loop
         Non_Zero_Printed := False;
         for J in S'Range(1) loop
            if J < 4 then
               -- rotation part:
               if S (I,J) /= 0.0 then
                  if S (I,J) < 0.0 then
                     Buffer (Pos) := '-';
                     Pos := Pos + 1;
                  else
                     if Non_Zero_Printed then
                        Buffer (Pos) := '+';
                        Pos := Pos + 1;
                     end if;
                  end if;
                  if Abs (S (I,J)) /= 1.0 then
                     declare
                        R : String := Rational_Translation (abs(S (I,J)));
                     begin
                        for Ch of R (2..R'Last) loop
                           Buffer (Pos) := Ch;
                           Pos := Pos + 1;
                        end loop;
                     end;
                     Buffer (Pos) := '*';
                     Pos := Pos + 1;
                  end if;
                  case J is
                     when 1 => Buffer (Pos) := 'X';
                     when 2 => Buffer (Pos) := 'Y';
                     when 3 => Buffer (Pos) := 'Z';
                     when others =>
                        raise CONSTRAINT_ERROR;
                  end case;
                  Pos := Pos + 1;
                  Non_Zero_Printed := True;
               end if;
            else
               -- translation part:
               for C of Rational_Translation (S (I,J)) loop
                  Buffer (Pos) := C;
                  Pos := Pos + 1;
               end loop;
            end if;
         end loop;
         if I < 3 then
            Buffer (Pos) := ',';
            Pos := Pos + 1;
         end if;
      end loop;
      
      return Buffer (1..Pos-1);
   end As_String;
   
begin
   
   if Exists ("DECODE_HALL_DEBUG") and then
     (
      Value ("DECODE_HALL_DEBUG") = "1" or else
        Value ("DECODE_HALL_DEBUG") = "true"
     )
   then
      Debug_Print_Matrices := True;
   end if;
   
   if Argument_Count > 0 and then 
     Index ("--help", Argument (1)) = 1
   then
      Put_Line ("Decode a Hall Crystallographic spacegroup symbol " &
                  "given on the command line");
      Put_Line ("and output symmtetry operators as general position " &
                  "point coordinates, e.g. '-X,-Y,Z+1/2");
      New_Line;
      Put_Line ("USAGE:");
      Put_Line ("  " & Command_Name & " 'P -2c'");
      Put_Line ("  " & Command_Name & " --help");
   else
      for I in 1 .. Argument_Count loop
         if Debug_Print_Matrices then
            Put_Line (Standard_Error, Argument (I));
         end if;
         
         declare
            Symmetry_Operators : Symmetry_Operator_Array :=
              Decode_Hall (Argument (I));
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
      end loop;
   end if;
   
end Decode_Hall;
