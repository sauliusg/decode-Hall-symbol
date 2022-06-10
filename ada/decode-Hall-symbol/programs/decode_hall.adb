with Text_IO;         use Text_IO;
with Ada.Integer_Text_IO; use Ada.Integer_Text_IO;
with Ada.Command_Line;    use Ada.Command_Line;
with Ada.Environment_Variables; use Ada.Environment_Variables;

procedure Decode_Hall is
   
   Debug_Print_Matrices : Boolean := False;
   
   UNKNOWN_AXIS : exception;
   UNKNOWN_ROTATION : exception;
   UNKNOWN_TRANSLATION : exception;
   
   type Axis_Direction_Type is (X_AXIS, Y_AXIS, Z_AXIS);
   
   type Axis_Order_Type is (IDENTITY, TWOFOLD, THREEFOLD, FOURFOLD, SIXFOLD);
   
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
      (4 =>  1.0, others => 0.0)
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
   
   procedure Put (F : File_Type; S : Symop) is
   begin
      for J in Symop'Range(1) loop
         for K in Symop'Range(1) loop
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
   
   Principal_Rotations : constant array 
     (Axis_Direction_Type, Axis_Order_Type) of Symop :=
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
   
   procedure Snap_To_Crystallographic_Translations (M : in out Symop) is
      Eps : constant Float := 1.0E-5;
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
         end if;
      end loop;
   end;
   
   procedure Add (M1 : in out Symop; M2 : in Symop) is
   begin
      for I in M1'Range(1) loop
         for J in M1'Range(2) loop
            M1 (I,J) := M1 (I,J) + M2(I,J);
         end loop;
      end loop;
      Snap_To_Crystallographic_Translations (M1);
   end;
   
   function "+" (M1, M2 : Symop) return Symop is
      M : Symop;
   begin
      for I in M1'Range(1) loop
         for J in M1'Range(2) loop
            M (I,J) := M1 (I,J) + M2(I,J);
         end loop;
      end loop;
      Snap_To_Crystallographic_Translations (M);
      return M;
   end;
   
   function "*" (M1, M2 : Symop) return Symop is
      M : Symop;
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
   end;
   
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
   
   function Rotation_Axis_Index (Rotation_Character : Character)
                                return Axis_Order_Type
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
   end;
   
   procedure Construct_Rotation_Matrix
     (
      Matrix : out Symop;
      Inversion : in Character;
      Axis : in Character;
      Rotation : in Character;
      Translations : String;
      Preceeding_Axis_Direction : in out Natural;
      Preceeding_Axis_Order : in out Natural;
      Axis_Number_Parameter : in Natural
     )
   is
      Axis_Number : Integer range 0..4 := Axis_Number_Parameter;
        
      function Ord (C : Character) return Positive is
         (Character'Pos (C) - Character'Pos ('0'));
      
   begin
      case Axis is 
         when 'x' =>
            Matrix :=
              Principal_Rotations (X_AXIS, Rotation_Axis_Index (Rotation));
            Axis_Number := 1;
         when 'y' =>
            Matrix :=
              Principal_Rotations (Y_AXIS, Rotation_Axis_Index (Rotation));
            Axis_Number := 2;
         when 'z' =>
            Matrix :=
              Principal_Rotations (Z_AXIS, Rotation_Axis_Index (Rotation));
            Axis_Number := 3;
         when ''' =>
            Matrix :=
              Face_Diagonal_Rotations (Preceeding_Axis_Direction, 1);
         when '"' =>
            Matrix :=
              Face_Diagonal_Rotations (Preceeding_Axis_Direction, 2);
         when '*' =>
            Matrix :=
              Body_Diagonal_Rotation;
         when others =>
            raise UNKNOWN_AXIS with "axis character " & Axis'Image;
      end case;
      
      Preceeding_Axis_Direction := Axis_Number;
      Preceeding_Axis_Order := Ord (Rotation);
      
      if Inversion = '-' then
         Matrix := Ci_Matrix * Matrix;
      end if;
      
      for Tr of Translations loop
         case Tr is
            when ' ' => null;
            when 'a' => Add (Matrix, To_Symop (Translation_a));
            when 'b' => Add (Matrix, To_Symop (Translation_b));
            when 'c' => Add (Matrix, To_Symop (Translation_c));
            when 'd' => Add (Matrix, To_Symop (Translation_d));
            when 'u' => Add (Matrix, To_Symop (Translation_u));
            when 'v' => Add (Matrix, To_Symop (Translation_v));
            when 'w' => Add (Matrix, To_Symop (Translation_w));
            when 'n' => Add (Matrix, To_Symop (Translation_n));
            when '1' => 
               case Rotation is 
                  when '3' => 
                     Add (Matrix, To_Symop (Translations_3_1, Axis_Number));
                  when '4' => 
                     Add (Matrix, To_Symop (Translations_4_1, Axis_Number));
                  when '6' => 
                     Add (Matrix, To_Symop (Translations_6_1, Axis_Number));
                  when others =>
                     raise UNKNOWN_ROTATION
                       with "mismatching translation " & 
                       Tr'Image & " for rotation " & Rotation'Image;
               end case;
            when '2' => 
               case Rotation is 
                  when '3' => 
                     Add (Matrix, To_Symop (Translations_3_2, Axis_Number));
                  when '6' => 
                     Add (Matrix, To_Symop (Translations_6_2, Axis_Number));
                  when others =>
                     raise UNKNOWN_ROTATION
                       with "mismatching translation " & 
                       Tr'Image & " for rotation " & Rotation'Image;
               end case;
            when '3' => 
               case Rotation is 
                  when '4' => 
                     Add (Matrix, To_Symop (Translations_4_3, Axis_Number));
                  when others =>
                     raise UNKNOWN_ROTATION
                       with "mismatching translation " & 
                       Tr'Image & " for rotation " & Rotation'Image;
               end case;
            when '4' => 
               case Rotation is 
                  when '6' => 
                     Add (Matrix, To_Symop (Translations_6_4, Axis_Number));
                  when others =>
                     raise UNKNOWN_ROTATION
                       with "mismatching translation " & 
                       Tr'Image & " for rotation " & Rotation'Image;
               end case;
            when '5' => 
               case Rotation is 
                  when '6' => 
                     Add (Matrix, To_Symop (Translations_6_5, Axis_Number));
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
      
   end;
      
   procedure Get_Hall_Symbol_Rotations
     (
      Symbol : in String;
      Pos : in out Positive;
      Rotations : out Symop_Array;
      N_Rotations : in out Natural;
      Preceeding_Axis_Direction : in out Natural;
      Preceeding_Axis_Order : in out Natural;
      Axis_Number : in Positive
     )
   is
      
      procedure Get_Inversion_Character (Inversion : out Character) is
      begin
         if Pos <= Symbol'Last and then
           Symbol (Pos) = '-' then
            Inversion := Symbol (Pos);
            Pos := Pos + 1;
         else
            Inversion := ' ';
         end if;
      end;
      
      procedure Get_Rotation_Character (Rotation : out Character) is
      begin
         if Pos <= Symbol'Last then
            pragma Assert (Symbol (Pos) in '2'..'6');
            Rotation := Symbol (Pos);
            Pos := Pos + 1;
         else
            Rotation := ' ';
         end if;
      end;
      
      procedure Get_Axis_Character (
                                    Axis : out Character;
                                    Axis_Number : in Positive;
                                    Rotation_Character : in Character;
                                    Preceeding_Axis_Order : in Natural
                                   ) is
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
                        if Preceeding_Axis_Order = 2 or else
                          Preceeding_Axis_Order = 4 then
                           Axis := 'x'; -- the \vec{a} crystallographic axis
                        elsif Preceeding_Axis_Order = 3 or else 
                          Preceeding_Axis_Order = 6 then
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
      end;
      
      procedure Get_Translation_Characters
        (
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
      end;      
      
      Inversion : Character;
      Rotation : Character;
      Axis : Character;
      Translations : String (1..3) := (others => ' ');
      N_Translations : Natural := 0;
      
   begin
      Skip_Spaces (Symbol, Pos);
      Get_Inversion_Character (Inversion);
      Get_Rotation_Character (Rotation);
      Get_Translation_Characters (Translations, N_Translations);
      Get_Axis_Character (Axis, Axis_Number, Rotation, Preceeding_Axis_Order);
      Get_Translation_Characters (Translations, N_Translations);
      
      if Rotation /= ' ' and then Axis /= ' ' then
         N_Rotations := N_Rotations + 1;
         Construct_Rotation_Matrix (Rotations (N_Rotations), 
                                    Inversion, Axis, Rotation,
                                    Translations, Preceeding_Axis_Direction,
                                    Preceeding_Axis_Order,
                                    Axis_Number);
         
         if Rotations (N_Rotations) = Unity_Matrix then
            N_Rotations := N_Rotations - 1;
         end if;
      end if;
   end Get_Hall_Symbol_Rotations;
   
   function Decode_Hall (Symbol : in String) return Symop_Array is
      Max_Symops : constant Integer := 192;
      
      Symops : Symop_Array (1 .. Max_Symops);
      N_Symops : Positive := 1;
      
      Pos : Positive := 1;      -- current position in the string 'Symbol'.
      
      N_Inversions : Positive;
      
      Centering : Symop_Array (1..4);
      N_Centering : Positive;
      
      Preceeding_Axis_Direction : Natural := 0;
      Preceeding_Axis_Order : Natural := 0;
      
      function Has_Symop
        (
         Symops : Symop_Array;
         Last_Symop_Index : Positive;
         Lookup_Symop : Symop
        )
        return Boolean is
      begin
         for I in 1 .. Last_Symop_Index loop
            if Symops (I) = Lookup_Symop then
               return True;
            end if;
         end loop;
         return False;
      end;
      
   begin
      Symops (1) := Unity_Matrix;
      
      Get_Hall_Symbol_Inversions (Symbol, Pos, N_Inversions);
      Get_Hall_Symbol_Centerings (Symbol, Pos, Centering, N_Centering);
      
      Get_Hall_Symbol_Rotations  (Symbol, Pos, Symops, N_Symops,
                                  Preceeding_Axis_Direction,
                                  Preceeding_Axis_Order, 1);

      Get_Hall_Symbol_Rotations  (Symbol, Pos, Symops, N_Symops, 
                                  Preceeding_Axis_Direction,
                                  Preceeding_Axis_Order, 2);

      Get_Hall_Symbol_Rotations  (Symbol, Pos, Symops, N_Symops,
                                  Preceeding_Axis_Direction,
                                  Preceeding_Axis_Order, 3);
      
      Get_Hall_Symbol_Rotations  (Symbol, Pos, Symops, N_Symops,
                                  Preceeding_Axis_Direction,
                                  Preceeding_Axis_Order, 4);
      
      if Debug_Print_Matrices then
         Put_Line (Standard_Error, "Inversions:");
         for I in 1..N_Inversions loop
            Put (Standard_Error, Inversion_Matrices (I));
            New_Line (Standard_Error);
         end loop;
         
         Put_Line (Standard_Error, "Centerings:");
         for I in 1..N_Centering loop
            Put (Standard_Error, Centering (I));
            New_Line (Standard_Error);
         end loop;
         
         Put_Line (Standard_Error, "Rotations:");
         for I in 1..N_Symops loop
            Put_Line (Standard_Error, I'Image);
            Put (Standard_Error, Symops (I));
            New_Line (Standard_Error);
         end loop;
      end if;
      
      -- Reconstruct all rotation operators:
      
      declare
         N, M : Positive := N_Symops;
         New_Symop : Symop;
      begin
         loop
            for I in 2 .. N loop
               New_Symop := Symops (I) * Symops (N);
               if not Has_Symop (Symops, M, New_Symop) then
                  M := M + 1;
                  Symops (M) := New_Symop;
               end if;
            end loop;
            N := N + 1;
            exit when N > M or else M >= Max_Symops;
         end loop;
         N_Symops := M;
      end;
      
      -- Add centering and inversion matrices:
      
      declare
         M : Positive := N_Symops;
         New_Symop : Symop;
      begin
         for I in 1..N_Inversions loop
            for C in 1..N_Centering loop
               if I /= 1 or else C /= 1 then
                  for S in 1..N_Symops loop
                     New_Symop :=
                       (Symops (S) + Centering (C)) * Inversion_Matrices (I);
                     if not Has_Symop (Symops, M, New_Symop) then
                       M := M + 1;
                       Symops (M) := New_Symop;
                     end if;
                  end loop;
               end if;
            end loop;
         end loop;
         N_Symops := M;
      end;      
      
      return Symops (1..N_Symops);
   end;
   
   function As_String (S : Symop) return String is
      Buffer : String (1..100); -- large enough to hold any symop
      Pos : Positive := 1;
      Non_Zero_Printed : Boolean;
      
      function Rational_Translation (T : Float) return String is
      begin
         if T = 0.0 then
            return "";
         elsif T = 0.5 then
            return "+1/2";
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
         else
            return T'Image;
         end if;
      end;
      
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
   end;
   
begin
   
   if Exists ("DECODE_HALL_DEBUG") and then
     (
      Value ("DECODE_HALL_DEBUG") = "1" or else
        Value ("DECODE_HALL_DEBUG") = "true"
     )
   then
      Debug_Print_Matrices := True;
   end if;
      
   for I in 1 .. Argument_Count loop
      if Debug_Print_Matrices then
         Put_Line (Standard_Error, Argument (I));
      end if;
      
      declare
         Symops : Symop_Array := Decode_Hall (Argument (I));
      begin
         if Debug_Print_Matrices then
            Put_Line (Standard_Error, "Symops:");
            for I in Symops'Range loop
               Put (Standard_Error, Symops (I));
               New_Line (Standard_Error);
            end loop;
         end if;
         
         for I in Symops'Range loop
            Put_Line (As_String (Symops (I)));
         end loop;
      end;
   end loop;
   
end Decode_Hall;
