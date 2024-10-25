with Ada.Text_IO; use Ada.Text_IO;

package Symmetry_Operations is
   
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
   
   -- Centering translations vectors from Hall 1981 [1], Table 1:
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
   
   -- Translation symbol vectors from Hall 1981 [1], Table 2, left side:
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
   
   -- Translation symbol vectors from Hall 1981 [1], Table 2, right
   --  side. These translations will have to be aplied along the
   --  specified axis:
   Translations_3_1 : constant Crystallographic_Translation_Component := (1,3);
   Translations_3_2 : constant Crystallographic_Translation_Component := (2,3);
   
   Translations_4_1 : constant Crystallographic_Translation_Component := (1,4);
   Translations_4_3 : constant Crystallographic_Translation_Component := (3,4);
   
   Translations_6_1 : constant Crystallographic_Translation_Component := (1,6);
   Translations_6_2 : constant Crystallographic_Translation_Component := (2,6);
   Translations_6_4 : constant Crystallographic_Translation_Component := (4,6);
   Translations_6_5 : constant Crystallographic_Translation_Component := (5,6);
   
   procedure Put (F : File_Type; S : Symmetry_Operator);
   
   procedure Put (F : File_Type; S : Crystallographic_Translation);
   
   function To_Symmetry_Operator (T : Crystallographic_Translation) 
                                 return Symmetry_Operator;
   
   function Axis_Index (Direction : Known_Axis_Direction) return Positive;
   
   function To_Symmetry_Operator (T : Crystallographic_Translation_Component;
                                  Axis_Direction : Known_Axis_Direction)
                                 return Symmetry_Operator;
   
   procedure Skip_Spaces (S : in String; Pos : in out Integer );

   procedure Decode_Centering_Symbol
     (
      Symbol : in String;
      Pos : in out Positive;
      Centering : out Symmetry_Operator_Array;
      N_Centering : out Positive
     );

end Symmetry_Operations;
