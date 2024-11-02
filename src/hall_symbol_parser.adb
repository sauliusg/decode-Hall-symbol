with Ada.Text_IO; use Ada.Text_IO;

with Parser_Tools;              use Parser_Tools;

with Symmetry_Operations;       use Symmetry_Operations;
with Change_Of_Basis;           use Change_Of_Basis;

package body Hall_Symbol_Parser is
   
   function IDENTITY return Axis_Order_Type 
     renames Symmetry_Operations.IDENTITY;
   
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
   -- 3. International Tables Volume B (2010), "Symmetry in reciprocal
   --  space". Section 1.4., Appendix A1.4.2. Space-group symbols for
   --  numeric and symbolic computations, URL:
   --  https://onlinelibrary.wiley.com/iucr/itc/B/ [accessed:
   --  2022-06-14T15:35+03:00]
   
   -- Rotation matrices from Hall 1981 [1], Table 3:
   
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
      
      Z_AXIS =>
        (-- axis z (c)
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
   
   -- Rotation matrices from Hall 1981 [1], Table 4:
   Face_Diagonal_Rotations :
     constant array (Known_Axis_Direction,1..2) of Symmetry_Operator :=
     (
      X_AXIS => (
                 1 => ( -- 2'
                       (-1.0,  0.0,  0.0,  0.0),
                       ( 0.0,  0.0, -1.0,  0.0),
                       ( 0.0, -1.0,  0.0,  0.0),
                       ( 0.0,  0.0,  0.0,  1.0)
                      ),
                 2 => ( -- 2"
                       (-1.0,  0.0,  0.0,  0.0),
                       ( 0.0,  0.0,  1.0,  0.0),
                       ( 0.0,  1.0,  0.0,  0.0),
                       ( 0.0,  0.0,  0.0,  1.0)
                      )
                ),
      Y_AXIS => (
                 1 => ( -- 2'
                       ( 0.0,  0.0, -1.0,  0.0),
                       ( 0.0, -1.0,  0.0,  0.0),
                       (-1.0,  0.0,  0.0,  0.0),
                       ( 0.0,  0.0,  0.0,  1.0)
                      ),
                 2 => ( -- 2"
                       ( 0.0,  0.0,  1.0,  0.0),
                       ( 0.0, -1.0,  0.0,  0.0),
                       ( 1.0,  0.0,  0.0,  0.0),
                       ( 0.0,  0.0,  0.0,  1.0)
                      )
                ),
      Z_AXIS => (
                 1 => ( -- 2'
                       ( 0.0, -1.0,  0.0,  0.0),
                       (-1.0,  0.0,  0.0,  0.0),
                       ( 0.0,  0.0, -1.0,  0.0),
                       ( 0.0,  0.0,  0.0,  1.0)
                      ),
                 2 => ( -- 2"
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
   
   -- -------------------------------------------------------------------------
   
   procedure Skip_Spaces (S : in String; Pos : in out Integer )
     renames Parser_Tools.Skip_Spaces;
   
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
     renames Symmetry_Operations.Decode_Centering_Symbol;
 
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
              with "rotation " & Character'Image (Rotation_Character);
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
            raise UNKNOWN_AXIS with "axis character " & Character'Image (Axis);
      end case;
      
      Preceeding_Axis_Order := Current_Axis_Order;
   end Get_Rotation_Matrix_From_Axis_And_Rotation;
   
   procedure Add_Translations_To_The_Rotation_Matrix
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
                       Character'Image (Tr) & " for rotation " & 
                       Character'Image (Rotation);
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
                       Character'Image (Tr) & " for rotation " &
                       Character'Image (Rotation);
               end case;
            when '3' => 
               case Rotation is 
                  when '4' => 
                     Add (Matrix, Translations_4_3, Axis_Direction);
                  when others =>
                     raise UNKNOWN_ROTATION
                       with "mismatching translation " & 
                       Character'Image (Tr) & " for rotation " &
                       Character'Image (Rotation);
               end case;
            when '4' => 
               case Rotation is 
                  when '6' => 
                     Add (Matrix, Translations_6_4, Axis_Direction);
                  when others =>
                     raise UNKNOWN_ROTATION
                       with "mismatching translation " & 
                       Character'Image (Tr) & " for rotation " &
                       Character'Image (Rotation);
               end case;
            when '5' => 
               case Rotation is 
                  when '6' => 
                     Add (Matrix, Translations_6_5, Axis_Direction);
                  when others =>
                     raise UNKNOWN_ROTATION
                       with "mismatching translation " & 
                       Character'Image (Tr) & " for rotation " &
                       Character'Image (Rotation);
               end case;
            when others =>
               raise UNKNOWN_TRANSLATION
                 with "translation character " & Character'Image (Tr);
         end case;
      end loop;
   end Add_Translations_To_The_Rotation_Matrix;
   
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
      
      Add_Translations_To_The_Rotation_Matrix
        (
         Matrix, Rotation,
         Translations,
         Axis_Direction
        );
   end Construct_Rotation_Matrix;
      
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
                          Axis_Order_Type'Image (Preceeding_Axis_Order) &
                          " and current axis " & 
                          Character'Image (Rotation_Character);
                     end if;
                  when '3' => 
                     Axis := '*';
                  when others =>
                     raise UNKNOWN_ROTATION 
                       with "wrong rotation character " & 
                       Character'Image (Rotation_Character) &
                       " for axis number" & Positive'Image (Axis_Number);
               end case;
            when 4 =>
               Axis := 'x';
            when others =>
               raise UNKNOWN_AXIS 
                 with "axis number" & Positive'Image (Axis_Number);
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
   
   procedure Apply_Change_Of_Basis
     (
      Symmetry_Operators : in out Symmetry_Operator_Array;
      N_Symmetry_Operators : in out Natural;

      Centering : in out Symmetry_Operator_Array;
      N_Centering : in out Natural;

      Inversions : in out Symmetry_Operator_Array;
      N_Inversions : in out Natural;

      Change_Of_Basis : in Symmetry_Operator
     )
   is separate;
   
   ----------------------------------------------------------------------------

   function Decode_Hall_Symbol (Symbol : in String)
                               return Symmetry_Operator_Array is
      Max_Symmetry_Operators : constant Integer := 192;
      
      Symmetry_Operators :
        Symmetry_Operator_Array (1 .. Max_Symmetry_Operators);
      N_Symmetry_Operators : Positive := 1;
      
      Pos : Positive := 1; -- current position in the string 'Symbol'.
      
      Inversions : Symmetry_Operator_Array (1 .. 2) :=
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
      Apply_Change_Of_Basis
        (
         Symmetry_Operators, N_Symmetry_Operators,
         Centering, N_Centering,
         Inversions, N_Inversions,
         Change_Of_Basis
        );
      
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
            Put_Line (Standard_Error, Integer'Image (I));
            Put (Standard_Error, Symmetry_Operators (I));
            New_Line (Standard_Error);
         end loop;
         
         Put_Line (Standard_Error, "Change of basis:");
         Put (Standard_Error, Change_Of_Basis);
         New_Line (Standard_Error);
      end if;
      
      -- Reconstruct all symmetry operators:
      
      Build_Group (Symmetry_Operators, N_Symmetry_Operators);
      
      return Symmetry_Operators (1..N_Symmetry_Operators);
   end Decode_Hall_Symbol;
   
end Hall_Symbol_Parser;
