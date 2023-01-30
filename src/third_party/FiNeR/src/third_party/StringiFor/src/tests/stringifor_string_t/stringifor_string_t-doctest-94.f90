program volatile_doctest
use stringifor_string_t
 type(string) :: astring
 character(len=:), allocatable :: acharacter
 logical :: test_passed(2)
 astring = '  one '
 acharacter = 'three'
 test_passed(1) = ((acharacter==astring).eqv..false.)
 astring = 'the same '
 acharacter = 'the same '
 test_passed(2) = ((acharacter==astring).eqv..true.)
 print '(L1)', all(test_passed)
endprogram volatile_doctest