MODULE strings

  USE shared_data

  IMPLICIT NONE

  PRIVATE :: integer4_as_string, integer8_as_string

  INTERFACE integer_as_string
    MODULE PROCEDURE integer4_as_string, integer8_as_string
  END INTERFACE integer_as_string

CONTAINS

  SUBROUTINE integer4_as_string(int_in, string)

    INTEGER(i4), INTENT(IN) :: int_in
    CHARACTER(LEN=*), INTENT(OUT) :: string

    INTEGER :: n_nums
    CHARACTER(LEN=9) :: numfmt

    IF (int_in == 0) THEN
      n_nums = 1
    ELSE
      n_nums = 1 + INT(LOG10(REAL(ABS(int_in), num)))
    ENDIF
    WRITE(numfmt, '(''(I'', I6.6, '')'')') n_nums
    WRITE(string, numfmt) int_in

  END SUBROUTINE integer4_as_string



  SUBROUTINE integer8_as_string(int_in, string)

    INTEGER(i8), INTENT(IN) :: int_in
    CHARACTER(LEN=*), INTENT(OUT) :: string

    INTEGER :: n_nums
    CHARACTER(LEN=12) :: numfmt

    IF (int_in == 0) THEN
      n_nums = 1
    ELSE
      n_nums = 1 + INT(LOG10(REAL(ABS(int_in), num)))
    ENDIF
    WRITE(numfmt, '(''(I'', I9.9, '')'')') n_nums
    WRITE(string, numfmt) int_in

  END SUBROUTINE integer8_as_string



  FUNCTION str_cmp(str_in, str_test)

    CHARACTER(*), INTENT(IN) ::  str_in, str_test
    CHARACTER(LEN=c_max_string_length) :: str_trim
    INTEGER :: test_len, in_len
    LOGICAL :: str_cmp

    str_trim = TRIM(ADJUSTL(str_in))
    test_len = LEN(TRIM(str_test))
    in_len = LEN(TRIM(str_trim))

    IF (test_len > 0) THEN
      IF (IACHAR(str_test(test_len:test_len)) == 0) test_len = test_len - 1
    ENDIF
    IF (in_len > 0) THEN
      IF (IACHAR(str_trim(in_len:in_len)) == 0) in_len = in_len - 1
    ENDIF

    IF (test_len /= in_len) THEN
      str_cmp = .FALSE.
      RETURN
    ENDIF

    str_cmp = (str_trim(1:test_len) == str_test(1:test_len))

  END FUNCTION str_cmp

END MODULE strings
