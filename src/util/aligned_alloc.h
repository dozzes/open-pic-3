#pragma once

#include <cstdlib>
#include <new>
#include <limits>

#if defined(_MSC_VER) || defined(__MINGW32__)
#include <malloc.h> // Necessary for MSVC-specific alignment functions
#endif

template <typename T, std::size_t Alignment> struct AlignedAllocator
{
    using value_type = T;

    // --- Key fix for MSVC: rebind support ---
    template <class U> struct rebind
    {
        using other = AlignedAllocator<U, Alignment>;
    };
    // -----------------------------------------

    AlignedAllocator() noexcept = default;
    template <class U> AlignedAllocator(const AlignedAllocator<U, Alignment>&) noexcept {}

    T* allocate(std::size_t n)
    {
        if (n == 0)
            return nullptr;
        if (n > std::numeric_limits<std::size_t>::max() / sizeof(T))
            throw std::bad_array_new_length();

        // Use the larger of the requested Alignment and the natural alignment of T.
        std::size_t actual_alignment = (Alignment > alignof(T)) ? Alignment : alignof(T);

        void* ptr = nullptr;
#if defined(_MSC_VER) || defined(__MINGW32__)
        // Use the computed alignment.
        ptr = _aligned_malloc(n * sizeof(T), actual_alignment);
        if (!ptr)
            throw std::bad_alloc();
#else
        if (posix_memalign(&ptr, actual_alignment, n * sizeof(T)) != 0)
            throw std::bad_alloc();
#endif
        return static_cast<T*>(ptr);
    }

    void deallocate(T* p, std::size_t) noexcept
    {
        if (!p)
            return;
#if defined(_MSC_VER) || defined(__MINGW32__)
        _aligned_free(p);
#else
        free(p);
#endif
    }
};

template <class T, class U, std::size_t Alignment>
bool operator==(const AlignedAllocator<T, Alignment>&, const AlignedAllocator<U, Alignment>&)
{
    return true;
}

template <class T, class U, std::size_t Alignment>
bool operator!=(const AlignedAllocator<T, Alignment>&, const AlignedAllocator<U, Alignment>&)
{
    return false;
}
