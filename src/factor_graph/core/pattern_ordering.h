#ifndef __PATTERN_ORDERING_H__
#define __PATTERN_ORDERING_H__

// for each variable, we need to know its dimension and its index in the H matrix.
struct PatternOrdering
{
public:
    PatternOrdering(const int& col = -1, const int& row = -1, const int& idx = -1)
        : m_col(col), m_row(row), m_idx(idx) {}

    inline PatternOrdering &SetCol(int val)
    {
        m_col = val;
        return *this;
    }

    inline PatternOrdering &SetRow(int val)
    {
        m_row = val;
        return *this;
    }

    inline PatternOrdering &SetIdx(int val)
    {
        m_idx = val;
        return *this;
    }

    inline int Row() const { return m_row; }
    inline int Col() const { return m_col; }
    inline int Idx() const { return m_idx; }
private:
    int m_row;
    int m_col;
    int m_idx; // The index of this variable in the H matrix.
};



#endif // !__PATTERN_ORDERING_H__