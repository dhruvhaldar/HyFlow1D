#pragma once

#include <iostream>
#include <string_view>

namespace Color {
    // Global flag to enable/disable colors
    inline bool enabled = true;

    // Wrapper struct for color codes
    struct Code {
        std::string_view value;
        // Implicit conversion to string_view for backward compatibility
        operator std::string_view() const { return value; }
    };

    // Overload operator<< to respect the enabled flag
    inline std::ostream& operator<<(std::ostream& os, const Code& c) {
        if (enabled) {
            os << c.value;
        }
        return os;
    }

    // Color Definitions
    inline constexpr Code Reset   {"\033[0m"};
    inline constexpr Code Bold    {"\033[1m"};
    inline constexpr Code Red     {"\033[31m"};
    inline constexpr Code Green   {"\033[32m"};
    inline constexpr Code Yellow  {"\033[33m"};
    inline constexpr Code Blue    {"\033[34m"};
    inline constexpr Code Magenta {"\033[35m"};
    inline constexpr Code Cyan    {"\033[36m"};
    inline constexpr Code White   {"\033[37m"};

    // Bold Variants
    inline constexpr Code BoldRed     {"\033[1;31m"};
    inline constexpr Code BoldGreen   {"\033[1;32m"};
    inline constexpr Code BoldYellow  {"\033[1;33m"};
    inline constexpr Code BoldBlue    {"\033[1;34m"};
    inline constexpr Code BoldMagenta {"\033[1;35m"};
    inline constexpr Code BoldCyan    {"\033[1;36m"};
    inline constexpr Code BoldWhite   {"\033[1;37m"};

    // Cursor Controls
    inline constexpr Code CursorHide  {"\033[?25l"};
    inline constexpr Code CursorShow  {"\033[?25h"};
}
