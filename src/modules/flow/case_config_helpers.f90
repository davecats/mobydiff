module case_config_helpers
    implicit none

    private

    public :: strip_config_comment, parse_config_section, split_config_key_value
    public :: to_lower, clean_config_string, next_config_entry

contains

    ! Read the next `key = value` entry from an open config file, skipping blank
    ! and comment lines and tracking the current [section]. ok is .false. at end
    ! of file. line_no counts every physical line read (for diagnostics). This is
    ! the shared scan skeleton for the per-case config readers.
    subroutine next_config_entry(unit, section, key, value, line_no, ok)
        integer, intent(in) :: unit
        character(len=*), intent(inout) :: section
        character(len=*), intent(out) :: key, value
        integer, intent(inout) :: line_no
        logical, intent(out) :: ok

        integer :: stat
        character(len=512) :: line

        ok = .false.
        do
            read(unit, '(A)', iostat=stat) line
            if (stat /= 0) return
            line_no = line_no + 1
            call strip_config_comment(line)
            line = adjustl(line)
            if (len_trim(line) == 0) cycle
            if (line(1:1) == "[") then
                call parse_config_section(line, section)
                cycle
            end if
            call split_config_key_value(line, key, value)
            ok = .true.
            return
        end do
    end subroutine next_config_entry

    subroutine strip_config_comment(line)
        character(len=*), intent(inout) :: line
        integer :: semicolon, hash, cut

        semicolon = index(line, ";")
        hash = index(line, "#")
        cut = 0
        if (semicolon > 0) cut = semicolon
        if (hash > 0 .and. (cut == 0 .or. hash < cut)) cut = hash
        if (cut > 0) line(cut:) = ""
    end subroutine strip_config_comment

    subroutine parse_config_section(line, section)
        character(len=*), intent(in) :: line
        character(len=*), intent(out) :: section
        integer :: last

        last = index(line, "]")
        if (last > 2) section = to_lower(trim(line(2:last-1)))
    end subroutine parse_config_section

    subroutine split_config_key_value(line, key, value)
        character(len=*), intent(in) :: line
        character(len=*), intent(out) :: key, value
        integer :: eq

        key = ""
        value = ""
        eq = index(line, "=")
        if (eq > 0) then
            key = trim(adjustl(line(:eq-1)))
            value = trim(adjustl(line(eq+1:)))
        else
            key = trim(adjustl(line))
        end if
    end subroutine split_config_key_value

    function to_lower(text) result(out)
        character(len=*), intent(in) :: text
        character(len=len(text)) :: out
        integer :: i, c

        do i = 1, len(text)
            c = iachar(text(i:i))
            if (c >= iachar("A") .and. c <= iachar("Z")) then
                out(i:i) = achar(c + iachar("a") - iachar("A"))
            else
                out(i:i) = text(i:i)
            end if
        end do
    end function to_lower

    function clean_config_string(text) result(out)
        character(len=*), intent(in) :: text
        character(len=len_trim(text)) :: out
        character(len=len_trim(text)) :: tmp
        integer :: n

        tmp = trim(adjustl(text))
        n = len_trim(tmp)
        if (n >= 2) then
            if ((tmp(1:1) == '"' .and. tmp(n:n) == '"') .or. &
                (tmp(1:1) == "'" .and. tmp(n:n) == "'")) then
                out = tmp(2:n-1)
                return
            end if
        end if
        out = tmp
    end function clean_config_string

end module case_config_helpers
